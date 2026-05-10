"""metrics/architecture.py
Exon and intron length statistics per transcript.

Output (<metrics_dir>/architecture.tsv) columns:
    transcript_id          joinable to manifest.tsv (preserved verbatim)
    gene_id                normalised; joinable across metrics
    strand
    n_exons
    n_internal_exons       max(0, n_exons - 2)
    first_exon_length      transcript direction; nt
    last_exon_length       transcript direction; nt; equals first when n_exons == 1
    internal_exon_mean     mean length of internal exons (nt); NA if n_internal_exons == 0
    internal_exon_median   median length of internal exons (nt); NA if n_internal_exons == 0
    internal_exon_sd       sample SD of internal exon lengths (nt); NA if n_internal_exons < 2
    intron_mean            mean intron length (genomic nt); NA if n_exons < 2
    intron_median          median intron length (genomic nt); NA if n_exons < 2
    intron_sd              sample SD of intron lengths (genomic nt); NA if n_exons < 3

Notes
-----
* First/last exon use transcript reading direction (so 'first exon' is at
  the 5' end of the mature mRNA, i.e. the exon with the highest genomic
  start on '-' strand).
* Intron lengths are genomic gaps between adjacent exons:
  next_exon.start - this_exon.end - 1, computed in genomic order. The
  resulting set is direction-independent so set-statistics (mean / median
  / SD) are unaffected by strand.
* SD is sample SD (Bessel-corrected, n-1 denominator) — matches R's sd().
  Defined only when the underlying set has at least two values.
* Single-exon transcripts: first_exon_length == last_exon_length, no
  introns, no internal exons; all stat columns are NA.
* Region junction counts and codon-relative distances live in
  junctions.tsv; this plugin deliberately stays at the architecture
  level (lengths and counts only).
"""
from __future__ import annotations

import csv
import logging
import statistics
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import (
    open_or_build_db,
    build_transcript_index,
    lookup_transcript,
    read_manifest_transcripts,
    normalise_id,
)

log = logging.getLogger('metrics.architecture')

OUTPUT_FILENAME = 'architecture.tsv'

_NA = 'NA'

_HEADER = [
    'transcript_id', 'gene_id', 'strand',
    'n_exons', 'n_internal_exons',
    'first_exon_length', 'last_exon_length',
    'internal_exon_mean', 'internal_exon_median', 'internal_exon_sd',
    'intron_mean', 'intron_median', 'intron_sd',
]


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [paths.canonical_gff, paths.manifest_tsv]


def compute(paths, metric_config, output_path: Path) -> None:
    """Entry point called by the orchestrator."""
    db = open_or_build_db(paths.canonical_gff,
                          paths.extract_dir / 'canonical.gff.db')

    index = build_transcript_index(db)
    if not index:
        log.error("No transcript-like features found in canonical.gff "
                  "(no features with exon children). Cannot proceed.")
        # Empty TSV with header; orchestrator's 'output produced' check passes.
        _write_tsv(output_path, [])
        return

    transcript_ids = read_manifest_transcripts(paths.manifest_tsv)
    log.info(f"Computing architecture for {len(transcript_ids)} transcripts")

    rows = []
    n_missing, n_no_exons = 0, 0
    sample_missing = []
    for tx_id in transcript_ids:
        tx_feature = lookup_transcript(db, tx_id, index)
        if tx_feature is None:
            n_missing += 1
            if len(sample_missing) < 5:
                sample_missing.append(tx_id)
            continue
        row = _compute_for_transcript(db, tx_feature, manifest_tx_id=tx_id)
        if row is None:
            n_no_exons += 1
            continue
        rows.append(row)

    if n_missing:
        log.warning(f"{n_missing} transcript(s) in manifest not found in canonical.gff")
        if n_missing >= max(5, len(transcript_ids) // 2):
            sample_keys = list(index.keys())[:8]
            log.warning("  ID format mismatch suspected. Examples:")
            log.warning(f"    Manifest IDs (first 5): {sample_missing}")
            log.warning(f"    GFF index keys (first 8): {sample_keys}")
            log.warning("  Run with -v for the full feature-type list. "
                        "If the formats clearly differ in a way the normaliser "
                        "doesn't handle, file a bug or extend KNOWN_PREFIXES "
                        "in lib/gff.py.")
    if n_no_exons:
        log.warning(f"{n_no_exons} transcript(s) had no exon features")

    _write_tsv(output_path, rows)
    log.info(f"Wrote {len(rows)} rows to {output_path.name}")


# ---------------------------------------------------------------------------
# Per-transcript computation
# ---------------------------------------------------------------------------

def _compute_for_transcript(db, transcript, manifest_tx_id: str) -> Optional[dict]:
    """Compute the row for one transcript. Returns None if no exons."""
    strand = transcript.strand
    exons = list(db.children(transcript, featuretype='exon'))
    if not exons:
        return None

    # Sort in transcript order to identify first / last / internal exons.
    if strand == '-':
        sorted_exons = sorted(exons, key=lambda e: e.start, reverse=True)
    else:
        sorted_exons = sorted(exons, key=lambda e: e.start)

    exon_lengths = [e.end - e.start + 1 for e in sorted_exons]
    n_exons = len(exon_lengths)

    first_len = exon_lengths[0]
    last_len = exon_lengths[-1]
    internal_lens = exon_lengths[1:-1]  # empty for n_exons <= 2

    # Intron lengths: genomic gaps between adjacent exons. Sort by genomic
    # start ascending — the resulting set is the same as a transcript-order
    # walk, and set-stats are direction-independent.
    by_start = sorted(exons, key=lambda e: e.start)
    intron_lens = [
        by_start[i + 1].start - by_start[i].end - 1
        for i in range(len(by_start) - 1)
    ]

    # gene_id from attributes; fall back to Parent or empty.
    gene_id = ''
    for attr in ('gene_id', 'Parent'):
        if attr in transcript.attributes and transcript.attributes[attr]:
            gene_id = normalise_id(transcript.attributes[attr][0])
            break

    return {
        'transcript_id': manifest_tx_id,   # preserve manifest ID exactly
        'gene_id': gene_id,
        'strand': strand,
        'n_exons': n_exons,
        'n_internal_exons': len(internal_lens),
        'first_exon_length': first_len,
        'last_exon_length': last_len,
        'internal_exon_mean': _mean(internal_lens),
        'internal_exon_median': _median(internal_lens),
        'internal_exon_sd': _stdev(internal_lens),
        'intron_mean': _mean(intron_lens),
        'intron_median': _median(intron_lens),
        'intron_sd': _stdev(intron_lens),
    }


# ---------------------------------------------------------------------------
# Stats helpers — return None for the under-defined cases (written as NA)
# ---------------------------------------------------------------------------

def _mean(values: list[int]) -> Optional[float]:
    if not values:
        return None
    return statistics.fmean(values)


def _median(values: list[int]):
    """Median of a list of ints. May return int (odd n) or float (even n)."""
    if not values:
        return None
    return statistics.median(values)


def _stdev(values: list[int]) -> Optional[float]:
    """Sample SD (n-1 denominator); matches R's sd(). Requires n >= 2."""
    if len(values) < 2:
        return None
    return statistics.stdev(values)


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
