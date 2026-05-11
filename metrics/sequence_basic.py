"""metrics/sequence_basic.py
Length and base composition per transcript per region.

Output (<metrics_dir>/sequence_basic.tsv) — long format, one row per
(transcript, region):

    transcript_id     joinable to manifest.tsv (preserved verbatim)
    gene_id           normalised; joinable across metrics
    region            region name (mRNA, CDS, 5UTR, 3UTR, ...)
    length            nt
    gc_content        (G+C) / length;          NA if length == 0
    frac_A, frac_C, frac_G, frac_U, frac_other   per-base fractions over length;
                                                  all NA if length == 0
    gc_skew           (G-C) / (G+C);          NA if G+C == 0
    at_skew           (A-U) / (A+U);          NA if A+U == 0
    purine_ratio      (A+G) / length;         NA if length == 0
    amino_ratio       (A+C) / length;         NA if length == 0
                                              (amino bases: A has NH2 at C6, C has NH2 at C4;
                                               keto bases G and U/T have C=O at the analogous
                                               position. IUPAC code M = aMino.)

Format
------
Long. Use pivot_wider() in R if region-specific columns are wanted
downstream. Filtering by `region` is the natural access pattern.

Notes
-----
* Region discovery: globs `extracted_*.fa` under extract_dir; region name
  is the filename's `<region>` middle. Uppercase canonicals are treated as
  ACGT/U; T and U both count as U (genome is typically DNA but the spec
  uses RNA convention). Anything else (N, R, Y, S, W, ...) is counted in
  `frac_other` and excluded from gc / skew numerators.
* UTR_pair format is a 3-line-per-record file with a ViennaRNA constraint
  line, not standard FASTA. Skipped by default (see _DEFAULT_SKIP).
* ID match between FASTA records and manifest: tries raw, prefix-stripped,
  fully normalised. Verbatim match should be the norm since extraction
  writes manifest IDs as headers, but the same defence as junctions.py
  is cheap.

Plugin config (all optional)
----------------------------
    regions: [<list>]         explicit region list; overrides auto-discovery.
                              Members map to extracted_<region>.fa.
    skip_regions: [<list>]    additional regions to skip beyond the default set.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import normalise_id, strip_prefix

log = logging.getLogger('metrics.sequence_basic')

OUTPUT_FILENAME = 'sequence_basic.tsv'

_NA = 'NA'

_HEADER = [
    'transcript_id', 'gene_id', 'region',
    'length', 'gc_content',
    'frac_A', 'frac_C', 'frac_G', 'frac_U', 'frac_other',
    'gc_skew', 'at_skew',
    'purine_ratio', 'amino_ratio',
]

# Regions whose extracted_*.fa files are not standard FASTA and must be
# skipped. UTR_pair is a 3-line format with a ViennaRNA constraint line.
_DEFAULT_SKIP: tuple[str, ...] = ('UTR_pair',)


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    inputs: list[Path] = [paths.manifest_tsv]
    inputs.extend(_discover_region_files(paths, metric_config).values())
    return inputs


def compute(paths, metric_config, output_path: Path) -> None:
    """Entry point called by the orchestrator."""
    region_files = _discover_region_files(paths, metric_config)
    if not region_files:
        log.error("No region FASTA files found under %s "
                  "(expected extracted_<region>.fa).", paths.extract_dir)
        _write_tsv(output_path, [])
        return

    # Manifest gives us the canonical transcript set + gene_id mapping.
    manifest = _read_manifest(paths.manifest_tsv)
    log.info(f"Manifest: {len(manifest)} transcripts. "
             f"Regions to process: {sorted(region_files)}")

    # Build a multi-key lookup from manifest IDs (raw / prefix-stripped /
    # normalised) → manifest row, so we can match FASTA headers tolerantly.
    manifest_index = _build_manifest_index(manifest)

    rows: list[dict] = []
    for region, fa_path in sorted(region_files.items()):
        log.info(f"Processing region '{region}' ({fa_path.name})")
        n_records, n_unmatched = 0, 0
        sample_unmatched: list[str] = []
        for record_id, sequence in _iter_fasta(fa_path):
            n_records += 1
            tx_row = _lookup_manifest(record_id, manifest_index)
            if tx_row is None:
                n_unmatched += 1
                if len(sample_unmatched) < 5:
                    sample_unmatched.append(record_id)
                continue
            rows.append(_compute_row(
                manifest_tx_id=tx_row['transcript_id'],
                gene_id=tx_row['gene_id'],
                region=region,
                sequence=sequence,
            ))
        log.info(f"  {n_records} records, {n_records - n_unmatched} matched manifest")
        if n_unmatched:
            log.warning(f"  {n_unmatched} FASTA record(s) not in manifest.")
            if n_unmatched >= max(5, n_records // 2):
                sample_keys = list(manifest_index.keys())[:8]
                log.warning("    Sample FASTA IDs: %s", sample_unmatched)
                log.warning("    Sample manifest keys: %s", sample_keys)

    _write_tsv(output_path, rows)
    log.info(f"Wrote {len(rows)} rows to {output_path.name}")


# ---------------------------------------------------------------------------
# Region / file discovery
# ---------------------------------------------------------------------------

def _discover_region_files(paths, metric_config) -> dict[str, Path]:
    """Resolve {region_name: fasta_path} for every region this plugin will read.

    If `metric_config['regions']` is provided, use it as the explicit list
    (still subject to skip_regions). Otherwise glob extract_dir.
    """
    extract_dir = paths.extract_dir
    skip = set(_DEFAULT_SKIP)
    skip.update(metric_config.get('skip_regions') or [])

    explicit = metric_config.get('regions')
    if explicit:
        files: dict[str, Path] = {}
        for region in explicit:
            if region in skip:
                continue
            fa = extract_dir / f'extracted_{region}.fa'
            files[region] = fa  # existence check happens at orchestrator level
        return files

    files = {}
    for fa in sorted(extract_dir.glob('extracted_*.fa')):
        # extracted_<region>.fa  ->  <region>
        region = fa.stem[len('extracted_'):]
        if not region or region in skip:
            continue
        files[region] = fa
    return files


# ---------------------------------------------------------------------------
# Manifest reading and lookup
# ---------------------------------------------------------------------------

def _read_manifest(manifest_tsv: Path) -> list[dict]:
    """Return list of {transcript_id, gene_id} dicts, one per manifest row.

    gene_id is normalised at read time; transcript_id is preserved verbatim.
    Missing gene_id column → empty string; missing values → empty string.
    """
    out = []
    with open(manifest_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'transcript_id' not in (reader.fieldnames or []):
            raise ValueError(
                f"manifest.tsv has no 'transcript_id' column; "
                f"got {reader.fieldnames}")
        has_gene = 'gene_id' in (reader.fieldnames or [])
        for row in reader:
            tid = row.get('transcript_id') or ''
            if not tid:
                continue
            gid_raw = row.get('gene_id', '') if has_gene else ''
            gid = normalise_id(gid_raw) if gid_raw else ''
            out.append({'transcript_id': tid, 'gene_id': gid})
    return out


def _build_manifest_index(manifest: list[dict]) -> dict[str, dict]:
    """Map every plausible form of a manifest transcript_id to its row.

    Three forms per ID: raw, prefix-stripped, fully normalised. First-write
    wins on collisions (shouldn't happen for a properly-filtered manifest).
    """
    index: dict[str, dict] = {}
    for row in manifest:
        tid = row['transcript_id']
        for key in (tid, strip_prefix(tid), normalise_id(tid)):
            if key:
                index.setdefault(key, row)
    return index


def _lookup_manifest(record_id: str, index: dict[str, dict]) -> Optional[dict]:
    for key in (record_id, strip_prefix(record_id), normalise_id(record_id)):
        if key in index:
            return index[key]
    return None


# ---------------------------------------------------------------------------
# FASTA reading — minimal parser, no biopython dependency required
# ---------------------------------------------------------------------------

def _iter_fasta(path: Path):
    """Yield (record_id, sequence_uppercase) for each record.

    record_id = first whitespace-delimited token of the header line.
    Sequences are uppercased and stripped of whitespace. Multi-line
    sequences are concatenated.
    """
    record_id: Optional[str] = None
    chunks: list[str] = []
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if record_id is not None:
                    yield record_id, ''.join(chunks)
                # Header: first whitespace-delimited token after '>'
                header = line[1:].strip()
                record_id = header.split(None, 1)[0] if header else ''
                chunks = []
            else:
                # Strip whitespace (incl. newline); uppercase
                chunks.append(line.strip().upper())
    if record_id is not None:
        yield record_id, ''.join(chunks)


# ---------------------------------------------------------------------------
# Per-record metric computation
# ---------------------------------------------------------------------------

def _compute_row(*, manifest_tx_id: str, gene_id: str, region: str,
                 sequence: str) -> dict:
    length = len(sequence)
    if length == 0:
        # Length-0 record: write 0 length, NA everywhere else.
        return {
            'transcript_id': manifest_tx_id,
            'gene_id': gene_id,
            'region': region,
            'length': 0,
            'gc_content': None,
            'frac_A': None, 'frac_C': None, 'frac_G': None,
            'frac_U': None, 'frac_other': None,
            'gc_skew': None, 'at_skew': None,
            'purine_ratio': None, 'amino_ratio': None,
        }

    # T and U are equivalent — bucket together. Anything else → 'other'.
    n_a = sequence.count('A')
    n_c = sequence.count('C')
    n_g = sequence.count('G')
    n_u = sequence.count('T') + sequence.count('U')
    n_other = length - (n_a + n_c + n_g + n_u)

    # GC content: (G+C) / length. Standard definition, includes all bases
    # in denominator so frac_other rolls into a transparent reduction of GC
    # without recomputation.
    gc_content = (n_g + n_c) / length

    # Skews — undefined when their respective denominators are zero.
    gc_skew = ((n_g - n_c) / (n_g + n_c)) if (n_g + n_c) > 0 else None
    at_skew = ((n_a - n_u) / (n_a + n_u)) if (n_a + n_u) > 0 else None

    # Directional ratios. Defined whenever length > 0 (handled above);
    # both have total length in the denominator, matching the frac_*
    # columns. Equal respectively to frac_A + frac_G and frac_A + frac_C
    # but exposed as named columns for analysis convenience.
    purine_ratio = (n_a + n_g) / length
    amino_ratio = (n_a + n_c) / length

    return {
        'transcript_id': manifest_tx_id,
        'gene_id': gene_id,
        'region': region,
        'length': length,
        'gc_content': gc_content,
        'frac_A': n_a / length,
        'frac_C': n_c / length,
        'frac_G': n_g / length,
        'frac_U': n_u / length,
        'frac_other': n_other / length,
        'gc_skew': gc_skew,
        'at_skew': at_skew,
        'purine_ratio': purine_ratio,
        'amino_ratio': amino_ratio,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _write_tsv(output_path: Path, rows: list[dict]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(
            f, fieldnames=_HEADER, delimiter='\t',
            extrasaction='ignore', lineterminator='\n',
        )
        writer.writeheader()
        for row in rows:
            out = {k: (_NA if row.get(k) is None else row.get(k)) for k in _HEADER}
            writer.writerow(out)
