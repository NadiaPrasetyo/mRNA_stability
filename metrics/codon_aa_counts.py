"""
codon_aa_counts — per-transcript codon and amino-acid composition.

Reads `extracted_CDS.fa`, walks codons in-frame from position 0, and
emits a wide-format table of raw counts: 64 standard codon columns
(`codon_AAA` … `codon_TTT`), a single `codon_other` column for codons
containing N or any non-ACGT base, 20 amino-acid columns (`aa_A` …
`aa_Y`), plus four derived summaries:

  - cds_length_codons : floor(len(CDS) / 3); total in-frame codons.
  - n_codons_scored   : codons counted in the 64 standard columns
                        (i.e. cds_length_codons - codon_other).
  - n_stops           : codon_TAA + codon_TAG + codon_TGA.
  - aa_total          : sum of the 20 aa_* columns
                        (i.e. n_codons_scored - n_stops).

Edge cases (agreed in handoff):
  - CDS length not a multiple of 3 → drop the incomplete final codon
    and log a warning.
  - Codon contains N or any non-ACGT base → bucketed into `codon_other`,
    not translated.
  - Stop codons counted in their own codon columns, excluded from all
    aa_* columns (terminal stop is present in codon counts, absent
    from AA totals).
  - U vs T: treated as equivalent at parse time.
  - Standard 64 codons and standard 20 AAs only (no Sec / Pyl).

FASTA records from `01_extract.py` carry composite headers of the form
`<gene_id>_<transcript_id>_<region>`. They are parsed with
`split_composite_fasta_id(..., region='CDS')` to recover the
transcript_id, then matched to the manifest via the same multi-form
lookup pattern (raw / prefix-stripped / fully normalised) used by
sequence_basic.py.

Transcripts present in `manifest.tsv` but absent from `extracted_CDS.fa`
(non-coding, or any other reason CDS extraction produced no record)
emit an NA row so the table preserves the manifest universe and joins
cleanly in R.

No reference-data dependency.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import (
    normalise_id,
    read_manifest_transcripts,
    split_composite_fasta_id,
    strip_prefix,
)

log = logging.getLogger('metrics.codon_aa_counts')

OUTPUT_FILENAME = 'codon_aa_counts.tsv'

_NA = 'NA'
_REGION = 'CDS'

# Standard genetic code, DNA / T form. Stops mapped to None.
_CODON_TABLE: dict[str, str | None] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": None,  "TAG": None,
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": None,  "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_STOP_CODONS: tuple[str, ...] = ('TAA', 'TAG', 'TGA')
_STANDARD_AAS: tuple[str, ...] = tuple('ACDEFGHIKLMNPQRSTVWY')    # 20, alpha
_ALL_CODONS: tuple[str, ...] = tuple(sorted(_CODON_TABLE.keys()))  # 64, alpha


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [
        paths.extract_dir / 'extracted_CDS.fa',
        paths.manifest_tsv,
    ]


def compute(paths, metric_config, output_path: Path) -> None:
    cds_fa = paths.extract_dir / 'extracted_CDS.fa'
    manifest_tsv = paths.manifest_tsv

    transcript_order = read_manifest_transcripts(manifest_tsv)
    manifest_index = _build_manifest_index(transcript_order)
    log.info(f"Manifest: {len(transcript_order)} transcripts.")

    rows_by_tid: dict[str, dict] = {}
    n_records = 0
    n_unmatched = 0
    sample_unmatched: list[str] = []

    for record_id, seq in _iter_fasta(cds_fa):
        n_records += 1
        manifest_tid = _resolve_to_manifest_tid(record_id, manifest_index)
        if manifest_tid is None:
            n_unmatched += 1
            if len(sample_unmatched) < 5:
                sample_unmatched.append(record_id)
            continue
        if manifest_tid in rows_by_tid:
            log.warning(
                "Duplicate match for manifest transcript %s "
                "(FASTA record %s); keeping first occurrence.",
                manifest_tid, record_id,
            )
            continue
        rows_by_tid[manifest_tid] = _count_one(seq, manifest_tid)

    if n_unmatched:
        log.warning(
            "%d/%d FASTA record(s) did not match manifest; sample: %s",
            n_unmatched, n_records, sample_unmatched,
        )

    header = (
        ['transcript_id']
        + [f'codon_{c}' for c in _ALL_CODONS]
        + ['codon_other']
        + [f'aa_{a}' for a in _STANDARD_AAS]
        + ['cds_length_codons', 'n_codons_scored', 'n_stops', 'aa_total']
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for tid in transcript_order:
            row = _format_row(tid, rows_by_tid.get(tid))
            out.write('\t'.join(row) + '\n')

    n_with_data = sum(1 for tid in transcript_order if tid in rows_by_tid)
    log.info(
        f"Wrote {len(transcript_order)} rows to {output_path.name} "
        f"({n_with_data} with data, {len(transcript_order) - n_with_data} NA)."
    )


# ---------------------------------------------------------------------------
# Manifest indexing and lookup (mirrors the pattern in sequence_basic.py;
# candidate for lifting to lib/gff.py once a third FASTA-consuming plugin
# lands).
# ---------------------------------------------------------------------------

def _build_manifest_index(transcript_ids: list[str]) -> dict[str, str]:
    """Map every plausible form of a manifest transcript_id to the canonical
    (verbatim) manifest tid. First-write wins on collisions.
    """
    index: dict[str, str] = {}
    for tid in transcript_ids:
        for key in (tid, strip_prefix(tid), normalise_id(tid)):
            if key:
                index.setdefault(key, tid)
    return index


def _resolve_to_manifest_tid(
    record_id: str,
    manifest_index: dict[str, str],
) -> Optional[str]:
    """Resolve a FASTA record header to the manifest's canonical
    transcript_id, returning None if no form matches.

    Stage 1: composite-header parse, anchored on the known region (CDS).
    Stage 2: try raw / prefix-stripped / fully-normalised forms of the
             parsed transcript portion, then of the original record_id
             as a final fallback (in case 01_extract.py was run without
             the composite-header convention).
    """
    parsed = split_composite_fasta_id(record_id, region=_REGION)
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
# FASTA reading — minimal parser, no biopython dependency
# ---------------------------------------------------------------------------

def _iter_fasta(fa_path: Path):
    """Yield (record_id, sequence) pairs from a FASTA file.

    record_id is the first whitespace-delimited token after `>`.
    Multi-line records and blank lines are handled. _count_one handles
    uppercase / U→T normalisation downstream.
    """
    header: Optional[str] = None
    chunks: list[str] = []
    with open(fa_path) as f:
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


# ---------------------------------------------------------------------------
# Per-transcript counts
# ---------------------------------------------------------------------------

def _count_one(seq: str, transcript_id: str) -> dict:
    """Walk codons in-frame from position 0; return per-transcript counts."""
    seq = seq.upper().replace('U', 'T')
    n = len(seq)
    trimmed = n - (n % 3)
    if n != trimmed:
        log.warning(
            "CDS length %d for %s is not a multiple of 3; dropping "
            "incomplete final codon (%d nt).",
            n, transcript_id, n - trimmed,
        )

    codon_counts: dict[str, int] = {c: 0 for c in _ALL_CODONS}
    aa_counts: dict[str, int] = {a: 0 for a in _STANDARD_AAS}
    codon_other = 0

    for i in range(0, trimmed, 3):
        codon = seq[i:i + 3]
        if codon in codon_counts:        # one of the 64 standard codons
            codon_counts[codon] += 1
            aa = _CODON_TABLE[codon]
            if aa is not None:           # exclude stops from AA tally
                aa_counts[aa] += 1
        else:
            codon_other += 1             # N / non-ACGT base in codon

    return {
        'codon_counts': codon_counts,
        'codon_other': codon_other,
        'aa_counts': aa_counts,
        'cds_length_codons': trimmed // 3,
    }


def _format_row(transcript_id: str, counts: Optional[dict]) -> list[str]:
    if counts is None:
        n_data_cols = len(_ALL_CODONS) + 1 + len(_STANDARD_AAS) + 4
        return [transcript_id] + [_NA] * n_data_cols

    codon_counts = counts['codon_counts']
    aa_counts = counts['aa_counts']
    codon_other = counts['codon_other']
    cds_length_codons = counts['cds_length_codons']

    n_codons_scored = sum(codon_counts.values())
    n_stops = sum(codon_counts[c] for c in _STOP_CODONS)
    aa_total = sum(aa_counts.values())

    row = [transcript_id]
    row += [str(codon_counts[c]) for c in _ALL_CODONS]
    row.append(str(codon_other))
    row += [str(aa_counts[a]) for a in _STANDARD_AAS]
    row += [
        str(cds_length_codons),
        str(n_codons_scored),
        str(n_stops),
        str(aa_total),
    ]
    return row