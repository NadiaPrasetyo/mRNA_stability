"""metrics/uorf.py
Per-transcript counts and summary statistics for upstream open reading
frames (uORFs) in the 5'UTR.

A uORF is defined here as a region beginning at an in-5'UTR start
codon (ATG by default; configurable) and ending at the first
downstream in-frame stop codon (TAA / TAG / TGA). Reading frame is
independent of the main ORF frame — uORFs in any of the three
possible frames are counted.

Two classes are distinguished:
  - classical uORF: stop codon falls entirely within the 5'UTR
  - overlapping uORF (oORF): stop codon falls at or past the main
    start codon, i.e., in the CDS or beyond

Output (<metrics_dir>/uorf.tsv) columns:
    transcript_id                       joinable to manifest.tsv
    gene_id                             normalised; joinable across metrics
    utr5_length                         nt; NA if no 5'UTR record
    n_uorfs                             count of classical uORFs
    n_overlapping_uorfs                 count of oORFs
    has_uorf                            0/1; 1 iff n_uorfs + n_overlapping_uorfs > 0
    total_classical_uorf_codons         sum of codon counts across classical uORFs;
                                            0 if none; NA if no 5'UTR record
    max_classical_uorf_codons           longest classical uORF (codons); NA if none
    dist_cap_to_first_uatg              nt from 5' end to first uATG (classical or
                                            overlapping); NA if no uATGs
    dist_last_uorf_stop_to_main_atg     intercistronic nt between the last classical
                                            uORF's stop and the main ATG; NA if no
                                            classical uORFs

Codon-count convention: includes the start codon, excludes the stop
codon (so "ATG TAA" = 1 codon, "ATG NNN TAA" = 2 codons).

Plugin config (all optional)
----------------------------
    start_codons: [<list>]   start codons to scan for; default ['ATG'].
                             Common extensions: CUG, GUG, ACG (near-cognate).
                             Both T and U forms accepted at parse time.
    min_codons: <int>        minimum codon count for a uORF to be counted;
                             default 1 (include single-codon ATG-stop).

Edge cases
----------
* Manifest transcript with no 5'UTR FASTA record (typically non-coding,
  or no annotated UTR): NA row, preserving manifest universe.
* 5'UTR present but shorter than 3 nt: zero counts (cannot host an ATG).
* uATG straddling the UTR-CDS boundary: not counted. uATG codon must
  fit entirely within the 5'UTR.
* uATG with no in-frame stop within (5'UTR + CDS): not counted.
  Logged at INFO level if it occurs. Could be extended to read into
  3'UTR if this matters in practice.
* uATGs in the same frame as the main ORF are counted normally — they
  produce either a classical uORF (if a stop occurs in the 5'UTR) or
  an oORF that typically terminates at the main stop codon.

No reference-data dependency.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import normalise_id, split_composite_fasta_id, strip_prefix

log = logging.getLogger('metrics.uorf')

OUTPUT_FILENAME = 'uorf.tsv'
_NA = 'NA'

_HEADER = [
    'transcript_id', 'gene_id',
    'utr5_length',
    'n_uorfs', 'n_overlapping_uorfs', 'has_uorf',
    'total_classical_uorf_codons', 'max_classical_uorf_codons',
    'dist_cap_to_first_uatg', 'dist_last_uorf_stop_to_main_atg',
]

_DEFAULT_START_CODONS: tuple[str, ...] = ('ATG',)
_DEFAULT_MIN_CODONS: int = 1
_STOP_CODONS: frozenset[str] = frozenset({'TAA', 'TAG', 'TGA'})


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [
        paths.manifest_tsv,
        paths.extract_dir / 'extracted_5UTR.fa',
        paths.extract_dir / 'extracted_CDS.fa',
    ]


def compute(paths, metric_config, output_path: Path) -> None:
    start_codons = _resolve_start_codons(metric_config)
    min_codons = int(metric_config.get('min_codons', _DEFAULT_MIN_CODONS))
    if min_codons < 1:
        log.warning(f"min_codons={min_codons} is below 1; clamping to 1.")
        min_codons = 1

    utr5_fa = paths.extract_dir / 'extracted_5UTR.fa'
    cds_fa = paths.extract_dir / 'extracted_CDS.fa'
    manifest_tsv = paths.manifest_tsv

    manifest = _read_manifest(manifest_tsv)
    manifest_index = _build_manifest_index(manifest)
    log.info(f"Manifest: {len(manifest)} transcripts. "
             f"Start codons: {sorted(start_codons)}. min_codons={min_codons}.")

    utr5_by_tid = _load_fasta_indexed(utr5_fa, manifest_index, region='5UTR')
    cds_by_tid = _load_fasta_indexed(cds_fa, manifest_index, region='CDS')

    rows: list[dict] = []
    n_no_utr = 0
    n_no_cds = 0
    n_unbounded = 0

    for row in manifest:
        tid = row['transcript_id']
        gid = row['gene_id']
        utr_seq = utr5_by_tid.get(tid)
        cds_seq = cds_by_tid.get(tid)

        if utr_seq is None:
            # Truly missing 5'UTR record: NA row to preserve manifest universe.
            n_no_utr += 1
            rows.append(_na_row(tid, gid))
            continue

        if cds_seq is None:
            # 5'UTR present but no CDS — pathological combination; emit NA.
            n_no_cds += 1
            rows.append(_na_row(tid, gid))
            continue

        stats = _compute_uorf_stats(
            utr_seq=utr_seq, cds_seq=cds_seq,
            start_codons=start_codons, min_codons=min_codons,
        )
        n_unbounded += stats.pop('n_unbounded')

        rows.append({'transcript_id': tid, 'gene_id': gid, **stats})

    if n_no_utr:
        log.info(f"{n_no_utr} manifest transcript(s) have no 5'UTR record "
                 f"(NA rows written).")
    if n_no_cds:
        log.warning(f"{n_no_cds} manifest transcript(s) have a 5'UTR record "
                    f"but no CDS record (NA rows written).")
    if n_unbounded:
        log.info(f"{n_unbounded} uATG(s) had no in-frame stop within the "
                 f"available 5'UTR+CDS sequence and were not counted.")

    _write_tsv(output_path, rows)
    log.info(f"Wrote {len(rows)} rows to {output_path.name}")


# ---------------------------------------------------------------------------
# Config resolution
# ---------------------------------------------------------------------------

def _resolve_start_codons(metric_config) -> frozenset[str]:
    raw = metric_config.get('start_codons') or list(_DEFAULT_START_CODONS)
    out = set()
    for c in raw:
        c = c.upper().replace('U', 'T')
        if len(c) != 3 or any(ch not in 'ACGT' for ch in c):
            log.warning(f"Skipping invalid start codon: {c!r}")
            continue
        out.add(c)
    if not out:
        log.warning("No valid start codons resolved; falling back to ATG.")
        return frozenset({'ATG'})
    return frozenset(out)


# ---------------------------------------------------------------------------
# Manifest reading and lookup (mirrors sequence_basic / codon_aa_counts)
# ---------------------------------------------------------------------------

def _read_manifest(manifest_tsv: Path) -> list[dict]:
    """Read manifest.tsv and return one row per transcript.

    The manifest is long-format with one row per (transcript, region).
    We filter to region == 'mRNA' to obtain the canonical
    transcript-level entry; every annotated transcript has an mRNA row.
    """
    out = []
    with open(manifest_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames or []
        for required in ('transcript_id', 'region'):
            if required not in fields:
                raise ValueError(
                    f"manifest.tsv missing required column {required!r}; "
                    f"got {fields}")
        has_gene = 'gene_id' in fields
        for row in reader:
            if row.get('region') != 'mRNA':
                continue
            tid = row.get('transcript_id') or ''
            if not tid:
                continue
            gid_raw = row.get('gene_id', '') if has_gene else ''
            gid = normalise_id(gid_raw) if gid_raw else ''
            out.append({'transcript_id': tid, 'gene_id': gid})
    return out


def _build_manifest_index(manifest: list[dict]) -> dict[str, str]:
    """Map every plausible form of a manifest transcript_id to the canonical
    (verbatim) manifest tid. First-write wins.
    """
    index: dict[str, str] = {}
    for row in manifest:
        tid = row['transcript_id']
        for key in (tid, strip_prefix(tid), normalise_id(tid)):
            if key:
                index.setdefault(key, tid)
    return index


def _resolve_to_manifest_tid(
    record_id: str,
    manifest_index: dict[str, str],
    region: str,
) -> Optional[str]:
    parsed = split_composite_fasta_id(record_id, region=region)
    candidates: list[str] = []
    if parsed is not None:
        candidates.append(parsed[1])
    candidates.append(record_id)

    for cand in candidates:
        for key in (cand, strip_prefix(cand), normalise_id(cand)):
            if key and key in manifest_index:
                return manifest_index[key]
    return None


# ---------------------------------------------------------------------------
# FASTA reading
# ---------------------------------------------------------------------------

def _iter_fasta(path: Path):
    header: Optional[str] = None
    chunks: list[str] = []
    with open(path) as f:
        for raw in f:
            line = raw.rstrip('\r\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)


def _load_fasta_indexed(
    path: Path,
    manifest_index: dict[str, str],
    region: str,
) -> dict[str, str]:
    """Load a FASTA into {manifest_tid: uppercase_T_form_sequence}.

    Records not matching any manifest transcript are dropped with a count.
    """
    out: dict[str, str] = {}
    n_records = 0
    n_unmatched = 0
    sample_unmatched: list[str] = []

    for record_id, seq in _iter_fasta(path):
        n_records += 1
        tid = _resolve_to_manifest_tid(record_id, manifest_index, region)
        if tid is None:
            n_unmatched += 1
            if len(sample_unmatched) < 5:
                sample_unmatched.append(record_id)
            continue
        seq = seq.upper().replace('U', 'T')
        if tid in out:
            log.warning(
                f"Duplicate match for manifest transcript {tid} "
                f"(FASTA record {record_id} in {path.name}); "
                f"keeping first occurrence.")
            continue
        out[tid] = seq

    if n_unmatched:
        log.warning(
            f"{n_unmatched}/{n_records} record(s) in {path.name} did not "
            f"match manifest; sample: {sample_unmatched}")

    return out


# ---------------------------------------------------------------------------
# uORF detection and summarisation
# ---------------------------------------------------------------------------

def _compute_uorf_stats(*, utr_seq: str, cds_seq: str,
                       start_codons: frozenset[str],
                       min_codons: int) -> dict:
    """Scan the 5'UTR for start codons; for each, find the first in-frame
    stop (continuing into the CDS if necessary); classify and summarise.

    Returns a dict containing all output columns plus 'n_unbounded'
    (uATGs with no in-frame stop found within 5'UTR + CDS). The caller
    pops 'n_unbounded' before writing the row.
    """
    utr_len = len(utr_seq)
    full = utr_seq + cds_seq
    full_len = len(full)

    n_uorfs = 0  # classical
    n_overlapping = 0
    n_unbounded = 0
    classical_codon_counts: list[int] = []
    first_uatg_pos: Optional[int] = None
    last_classical_stop_end: Optional[int] = None  # pos one past the last classical stop

    # Iterate every 5'UTR position that could host a complete start codon.
    # range(utr_len - 2) gives i ∈ [0, utr_len - 3], so the codon at
    # full[i:i+3] is guaranteed to lie entirely within the 5'UTR.
    for i in range(utr_len - 2):
        codon = full[i:i + 3]
        if codon not in start_codons:
            continue

        # Find first in-frame stop. range(i+3, full_len-2, 3) gives j
        # positions where j+3 <= full_len (codon fits in full sequence).
        stop_pos = None
        for j in range(i + 3, full_len - 2, 3):
            if full[j:j + 3] in _STOP_CODONS:
                stop_pos = j
                break

        if stop_pos is None:
            n_unbounded += 1
            continue

        codon_count = (stop_pos - i) // 3  # includes start, excludes stop
        if codon_count < min_codons:
            continue

        if first_uatg_pos is None or i < first_uatg_pos:
            first_uatg_pos = i

        if stop_pos + 3 <= utr_len:
            # Stop codon ends within the 5'UTR → classical uORF.
            n_uorfs += 1
            classical_codon_counts.append(codon_count)
            stop_end = stop_pos + 3
            if (last_classical_stop_end is None
                    or stop_end > last_classical_stop_end):
                last_classical_stop_end = stop_end
        else:
            # Stop codon extends to or past the main ATG → overlapping uORF.
            n_overlapping += 1

    total_classical = sum(classical_codon_counts)
    max_classical = max(classical_codon_counts) if classical_codon_counts else None
    dist_last_stop_to_main = (utr_len - last_classical_stop_end
                              if last_classical_stop_end is not None else None)
    has_uorf = 1 if (n_uorfs + n_overlapping) > 0 else 0

    return {
        'utr5_length': utr_len,
        'n_uorfs': n_uorfs,
        'n_overlapping_uorfs': n_overlapping,
        'has_uorf': has_uorf,
        'total_classical_uorf_codons': total_classical,
        'max_classical_uorf_codons': max_classical,
        'dist_cap_to_first_uatg': first_uatg_pos,
        'dist_last_uorf_stop_to_main_atg': dist_last_stop_to_main,
        'n_unbounded': n_unbounded,
    }


def _na_row(transcript_id: str, gene_id: str) -> dict:
    return {
        'transcript_id': transcript_id,
        'gene_id': gene_id,
        'utr5_length': None,
        'n_uorfs': None,
        'n_overlapping_uorfs': None,
        'has_uorf': None,
        'total_classical_uorf_codons': None,
        'max_classical_uorf_codons': None,
        'dist_cap_to_first_uatg': None,
        'dist_last_uorf_stop_to_main_atg': None,
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
