"""metrics/stopfree.py
Longest stop-codon-free segment per transcript region.

Scans each FASTA record independently and reports the longest nucleotide span
that contains no in-frame stop codon in the tested reading frames. By default
only the transcript-forward strand is scanned, because extracted transcript
regions are already written in 5'→3' transcript orientation. Set
``forward_only: false`` in the dataset YAML to scan all six frames.

Output (<metrics_dir>/stopfree.tsv) — long format, one row per
(transcript, region):

    transcript_id             joinable to manifest.tsv (preserved verbatim)
    gene_id                   normalised; joinable across metrics
    region                    region name (mRNA, CDS, 5UTR, 3UTR, ...)
    sequence_length           nt
    stopfree_length           nt; longest span without an in-frame stop codon
    stopfree_fraction         stopfree_length / sequence_length; NA if length 0
    direction                 forward or reverse
    frame                     1, 2, or 3 relative to the analysed orientation
    stopfree_start            1-based inclusive coordinate in analysed orientation
    stopfree_end              1-based inclusive coordinate in analysed orientation
    source_start              1-based inclusive coordinate in the input sequence
    source_end                1-based inclusive coordinate in the input sequence
    n_full_codons_scanned     number of full codon positions assessed
    n_stop_codons             number of stop codons found in assessed positions

Coordinates
-----------
For forward hits, stopfree_* and source_* coordinates are identical. For
reverse hits, stopfree_* coordinates refer to the reverse-complement sequence;
source_* maps the same interval back to the original input sequence.

A trailing partial codon is included in the final stop-free span. The scan only
uses complete codons to detect stops, but the metric itself is a nucleotide
span until the next in-frame stop codon rather than a complete-codon count.

Plugin config (all optional)
----------------------------
    regions: [<list>]         explicit region list; overrides auto-discovery.
                              Members map to extracted_<region>.fa.
    skip_regions: [<list>]    additional regions to skip beyond UTR_pair.
    forward_only: true        true = 3 transcript-forward frames;
                              false = all 6 frames.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import normalise_id, strip_prefix, split_composite_fasta_id

log = logging.getLogger('metrics.stopfree')

OUTPUT_FILENAME = 'stopfree.tsv'

_NA = 'NA'
_STOP_CODONS = {'TAA', 'TAG', 'TGA'}
_DEFAULT_SKIP: tuple[str, ...] = ('UTR_pair',)
_RC_TRANS = str.maketrans({
    # Canonical DNA/RNA bases. U is treated as T before complementing.
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
    'a': 'T', 'c': 'G', 'g': 'C', 't': 'A', 'u': 'A',
    # Ambiguous IUPAC bases are retained as unknown rather than raising.
    'N': 'N', 'R': 'N', 'Y': 'N', 'S': 'N', 'W': 'N', 'K': 'N',
    'M': 'N', 'B': 'N', 'D': 'N', 'H': 'N', 'V': 'N',
    'n': 'N', 'r': 'N', 'y': 'N', 's': 'N', 'w': 'N', 'k': 'N',
    'm': 'N', 'b': 'N', 'd': 'N', 'h': 'N', 'v': 'N',
})

_HEADER = [
    'transcript_id', 'gene_id', 'region',
    'sequence_length', 'stopfree_length', 'stopfree_fraction',
    'direction', 'frame',
    'stopfree_start', 'stopfree_end',
    'source_start', 'source_end',
    'n_full_codons_scanned', 'n_stop_codons',
]


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    inputs: list[Path] = [paths.manifest_tsv]
    inputs.extend(_discover_region_files(paths, metric_config).values())
    return inputs


def compute(paths, metric_config, output_path: Path) -> None:
    region_files = _discover_region_files(paths, metric_config)
    if not region_files:
        log.error(
            "No region FASTA files found under %s "
            "(expected extracted_<region>.fa).",
            paths.extract_dir,
        )
        _write_tsv(output_path, [])
        return

    forward_only = bool(metric_config.get('forward_only', True))
    manifest = _read_manifest(paths.manifest_tsv)
    manifest_index = _build_manifest_index(manifest)

    log.info(
        "Manifest: %d transcripts. Regions to process: %s. Mode: %s.",
        len(manifest), sorted(region_files),
        'forward-only' if forward_only else 'six-frame',
    )

    rows: list[dict] = []
    for region, fa_path in sorted(region_files.items()):
        log.info("Processing region '%s' (%s)", region, fa_path.name)
        n_records, n_unmatched = 0, 0
        sample_unmatched: list[str] = []

        for record_id, sequence in _iter_fasta(fa_path):
            n_records += 1
            tx_row = _lookup_manifest(record_id, manifest_index, region=region)
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
                forward_only=forward_only,
            ))

        log.info("  %d records, %d matched manifest", n_records, n_records - n_unmatched)
        if n_unmatched:
            log.warning("  %d FASTA record(s) not in manifest.", n_unmatched)
            if n_unmatched >= max(5, n_records // 2):
                log.warning("    Sample FASTA IDs: %s", sample_unmatched)
                log.warning("    Sample manifest keys: %s", list(manifest_index.keys())[:8])

    _write_tsv(output_path, rows)
    log.info("Wrote %d rows to %s", len(rows), output_path.name)


# ---------------------------------------------------------------------------
# Region / file discovery
# ---------------------------------------------------------------------------

def _discover_region_files(paths, metric_config) -> dict[str, Path]:
    extract_dir = paths.extract_dir
    skip = set(_DEFAULT_SKIP)
    skip.update(metric_config.get('skip_regions') or [])

    explicit = metric_config.get('regions')
    if explicit:
        files: dict[str, Path] = {}
        for region in explicit:
            if region in skip:
                continue
            files[region] = extract_dir / f'extracted_{region}.fa'
        return files

    files = {}
    for fa in sorted(extract_dir.glob('extracted_*.fa')):
        region = fa.stem[len('extracted_'):]
        if not region or region in skip:
            continue
        files[region] = fa
    return files


# ---------------------------------------------------------------------------
# Manifest reading and lookup
# ---------------------------------------------------------------------------

def _read_manifest(manifest_tsv: Path) -> list[dict]:
    """Read manifest.tsv and return one transcript-level row per transcript."""
    out = []
    with open(manifest_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        for required in ('transcript_id', 'region'):
            if required not in fields:
                raise ValueError(
                    f"manifest.tsv missing required column {required!r}; got {fields}"
                )
        has_gene = 'gene_id' in fields
        for row in reader:
            if row.get('region') != 'mRNA':
                continue
            tid = row.get('transcript_id') or ''
            if not tid:
                continue
            gid_raw = row.get('gene_id', '') if has_gene else ''
            out.append({
                'transcript_id': tid,
                'gene_id': normalise_id(gid_raw) if gid_raw else '',
            })
    return out


def _build_manifest_index(manifest: list[dict]) -> dict[str, dict]:
    index: dict[str, dict] = {}
    for row in manifest:
        tid = row['transcript_id']
        for key in (tid, strip_prefix(tid), normalise_id(tid)):
            if key:
                index.setdefault(key, row)
    return index


def _lookup_manifest(record_id: str, index: dict[str, dict], region: str) -> Optional[dict]:
    parsed = split_composite_fasta_id(record_id, region=region)
    candidates: list[str] = []
    if parsed is not None:
        candidates.append(parsed[1])
    candidates.extend((record_id, strip_prefix(record_id), normalise_id(record_id)))

    seen = set()
    for cand in candidates:
        if not cand or cand in seen:
            continue
        seen.add(cand)
        if cand in index:
            return index[cand]
    return None


# ---------------------------------------------------------------------------
# FASTA and sequence helpers
# ---------------------------------------------------------------------------

def _iter_fasta(path: Path):
    """Yield (record_id, sequence) for each FASTA record in a multi-FASTA file."""
    record_id: Optional[str] = None
    chunks: list[str] = []
    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if record_id is not None:
                    yield record_id, ''.join(chunks)
                header = line[1:].strip()
                record_id = header.split(None, 1)[0] if header else ''
                chunks = []
            else:
                chunks.append(line)
        if record_id is not None:
            yield record_id, ''.join(chunks)


def _clean_sequence(sequence: str) -> str:
    return sequence.upper().replace('U', 'T')


def reverse_complement(sequence: str) -> str:
    """Reverse-complement a DNA/RNA-ish sequence, preserving unknown bases as N."""
    return sequence.translate(_RC_TRANS)[::-1].upper().replace('U', 'T')


# ---------------------------------------------------------------------------
# Stop-free scan
# ---------------------------------------------------------------------------

def longest_stopfree_region(sequence: str, forward_only: bool = True) -> dict:
    """Return the longest stop-free interval across 3 or 6 reading frames.

    Internal coordinates are 0-based, half-open in the analysed orientation.
    Output formatting converts them to 1-based inclusive coordinates.
    """
    seq_forward = _clean_sequence(sequence)
    sequence_length = len(seq_forward)

    if sequence_length == 0:
        return {
            'sequence_length': 0,
            'direction': None,
            'frame0': None,
            'start0': None,
            'end0': None,
            'length': 0,
            'source_start0': None,
            'source_end0': None,
            'n_full_codons_scanned': 0,
            'n_stop_codons': 0,
        }

    directions = [('forward', seq_forward)]
    if not forward_only:
        directions.append(('reverse', reverse_complement(seq_forward)))

    best: Optional[dict] = None
    n_full_codons_scanned = 0
    n_stop_codons = 0

    for direction, seq in directions:
        for frame0 in range(3):
            region_start = min(frame0, sequence_length)

            # If the sequence is shorter than this frame offset, the frame has
            # an empty interval. Keeping it non-negative avoids the standalone
            # script's short-sequence oddity while preserving deterministic ties.
            if frame0 >= sequence_length:
                candidate = {
                    'direction': direction,
                    'frame0': frame0,
                    'start0': sequence_length,
                    'end0': sequence_length,
                    'length': 0,
                }
                best = _choose_better(best, candidate)
                continue

            for i in range(frame0, sequence_length - 2, 3):
                n_full_codons_scanned += 1
                codon = seq[i:i + 3]
                if codon in _STOP_CODONS:
                    n_stop_codons += 1
                    candidate = {
                        'direction': direction,
                        'frame0': frame0,
                        'start0': region_start,
                        'end0': i,
                        'length': i - region_start,
                    }
                    best = _choose_better(best, candidate)
                    region_start = i + 3

            candidate = {
                'direction': direction,
                'frame0': frame0,
                'start0': region_start,
                'end0': sequence_length,
                'length': sequence_length - region_start,
            }
            best = _choose_better(best, candidate)

    assert best is not None
    source_start0, source_end0 = _map_to_source_coords(
        start0=best['start0'],
        end0=best['end0'],
        sequence_length=sequence_length,
        direction=best['direction'],
    )

    return {
        'sequence_length': sequence_length,
        'direction': best['direction'],
        'frame0': best['frame0'],
        'start0': best['start0'],
        'end0': best['end0'],
        'length': best['length'],
        'source_start0': source_start0,
        'source_end0': source_end0,
        'n_full_codons_scanned': n_full_codons_scanned,
        'n_stop_codons': n_stop_codons,
    }


def _choose_better(best: Optional[dict], candidate: dict) -> dict:
    """Prefer the longest interval; stable ties keep first direction/frame/span."""
    if best is None:
        return candidate
    if candidate['length'] > best['length']:
        return candidate
    return best


def _map_to_source_coords(*, start0: int, end0: int,
                          sequence_length: int, direction: str) -> tuple[int, int]:
    """Map a half-open interval in analysed orientation to source coordinates.

    Returned coordinates remain 0-based, half-open. For reverse-complement
    intervals, the mapped source interval is normalised so start <= end.
    """
    if direction == 'forward':
        return start0, end0

    # Reverse complement coordinate r maps to source coordinate n - 1 - r.
    # A half-open interval [start0, end0) on RC spans source coordinates
    # [n - end0, n - start0) in forward orientation.
    return sequence_length - end0, sequence_length - start0


# ---------------------------------------------------------------------------
# Per-record row formatting
# ---------------------------------------------------------------------------

def _compute_row(*, manifest_tx_id: str, gene_id: str, region: str,
                 sequence: str, forward_only: bool) -> dict:
    result = longest_stopfree_region(sequence, forward_only=forward_only)
    seq_len = result['sequence_length']
    stopfree_len = result['length']

    return {
        'transcript_id': manifest_tx_id,
        'gene_id': gene_id,
        'region': region,
        'sequence_length': seq_len,
        'stopfree_length': stopfree_len,
        'stopfree_fraction': (stopfree_len / seq_len) if seq_len else None,
        'direction': result['direction'],
        'frame': (result['frame0'] + 1) if result['frame0'] is not None else None,
        'stopfree_start': _to_1based_start(result['start0'], stopfree_len),
        'stopfree_end': _to_1based_end(result['end0'], stopfree_len),
        'source_start': _to_1based_start(result['source_start0'], stopfree_len),
        'source_end': _to_1based_end(result['source_end0'], stopfree_len),
        'n_full_codons_scanned': result['n_full_codons_scanned'],
        'n_stop_codons': result['n_stop_codons'],
    }


def _to_1based_start(start0: Optional[int], length: int) -> Optional[int]:
    if start0 is None or length == 0:
        return None
    return start0 + 1


def _to_1based_end(end0: Optional[int], length: int) -> Optional[int]:
    if end0 is None or length == 0:
        return None
    return end0


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
