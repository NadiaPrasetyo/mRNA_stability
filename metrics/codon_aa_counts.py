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

Transcripts present in `manifest.tsv` but absent from `extracted_CDS.fa`
(non-coding, or any other reason CDS extraction produced no record)
emit an NA row so the table preserves the manifest universe and joins
cleanly in R.

No reference-data dependency.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

from lib.gff import read_manifest_transcripts

logger = logging.getLogger(__name__)

OUTPUT_FILENAME = "codon_aa_counts.tsv"

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

_STOP_CODONS: tuple[str, ...] = ("TAA", "TAG", "TGA")
_STANDARD_AAS: tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWY")    # 20, alpha
_ALL_CODONS: tuple[str, ...] = tuple(sorted(_CODON_TABLE.keys()))  # 64, alpha


def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [
        paths.extract_dir / "extracted_CDS.fa",
        paths.extract_dir / "manifest.tsv",
    ]


def _iter_fasta(fa_path: Path):
    """Yield (id, sequence) pairs from a FASTA file.

    ID is the first whitespace-delimited token after `>`. Multi-line
    records and blank lines are handled.
    """
    header: str | None = None
    chunks: list[str] = []
    with open(fa_path) as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                # Header: first whitespace-delimited token after '>'
                header = line[1:].strip()
                record_id = header.split(None, 1)[0] if header else ''
                
                # Patch: Extract transcript_id from <gene>_<transcript>_<region>
                if record_id and record_id.count('_') >= 2:
                    # Assuming format ENSG..._ENST..._REGION
                    parts = record_id.split('_')
                    record_id = parts[1]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def _count_one(seq: str, transcript_id: str) -> dict:
    """Walk codons in-frame from position 0; return per-transcript counts."""
    seq = seq.upper().replace("U", "T")
    n = len(seq)
    trimmed = n - (n % 3)
    if n != trimmed:
        logger.warning(
            "CDS length %d for %s is not a multiple of 3; dropping "
            "incomplete final codon (%d nt)",
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
        "codon_counts": codon_counts,
        "codon_other": codon_other,
        "aa_counts": aa_counts,
        "cds_length_codons": trimmed // 3,
    }


def _format_row(transcript_id: str, counts: dict | None) -> list[str]:
    if counts is None:
        # Manifest transcript with no CDS record — emit NA row to
        # preserve the manifest universe.
        n_data_cols = len(_ALL_CODONS) + 1 + len(_STANDARD_AAS) + 4
        return [transcript_id] + ["NA"] * n_data_cols

    codon_counts = counts["codon_counts"]
    aa_counts = counts["aa_counts"]
    codon_other = counts["codon_other"]
    cds_length_codons = counts["cds_length_codons"]

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


def compute(paths, metric_config, output_path: Path) -> None:
    cds_fa = paths.extract_dir / "extracted_CDS.fa"
    manifest_tsv = paths.extract_dir / "manifest.tsv"

    transcript_order = read_manifest_transcripts(manifest_tsv)

    rows_by_id: dict[str, dict] = {}
    for fasta_id, seq in _iter_fasta(cds_fa):
        if fasta_id in rows_by_id:
            logger.warning(
                "Duplicate transcript id %s in %s; keeping first occurrence",
                fasta_id, cds_fa.name,
            )
            continue
        rows_by_id[fasta_id] = _count_one(seq, fasta_id)

    extra = set(rows_by_id) - set(transcript_order)
    if extra:
        logger.warning(
            "%d transcript id(s) in %s absent from manifest (will not "
            "appear in output); example: %s",
            len(extra), cds_fa.name, next(iter(extra)),
        )

    header = (
        ["transcript_id"]
        + [f"codon_{c}" for c in _ALL_CODONS]
        + ["codon_other"]
        + [f"aa_{a}" for a in _STANDARD_AAS]
        + ["cds_length_codons", "n_codons_scored", "n_stops", "aa_total"]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for tid in transcript_order:
            row = _format_row(tid, rows_by_id.get(tid))
            out.write("\t".join(row) + "\n")
