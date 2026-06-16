"""lib/fasta_to_extract.py
Generate the extract_dir file tree that GFF-dependent metrics plugins
expect, sourced entirely from an annotated FASTA file.

For WT sequences the real pipeline populates extract_dir by running a
prior extraction step against a reference GFF.  For engineered (or
any FASTA-only) sequences, this module generates the same files from
header-embedded annotations so the plugins themselves need zero
modification.

Files produced under <extract_dir>/:
  manifest.tsv          ─ one mRNA row per sequence (+ 5UTR, CDS, 3UTR rows)
  canonical.gff         ─ GFF3 with exon/CDS/UTR features derived from headers
  extracted_5UTR.fa     ─ FASTA of 5' UTR sequences (absent if no 5' UTR)
  extracted_CDS.fa      ─ FASTA of CDS sequences
  extracted_3UTR.fa     ─ FASTA of 3' UTR sequences (absent if no 3' UTR)
  extracted_mRNA.fa     ─ full transcript sequences (copy of input FASTA)

The generated canonical.gff uses the sequence ID as both seqid and
transcript/gene ID so gffutils can build its DB normally.  Exon
features are synthesised from the EXONS= tag when present, or from the
CDS= tag alone (producing a single-exon transcript) when EXONS= is
absent.

Public API
──────────
  build_extract_dir(
      fasta_path : Path,
      extract_dir: Path,
      force      : bool = False,
  ) -> ExtractSummary

ExtractSummary.n_sequences    total sequences processed
ExtractSummary.n_with_cds     sequences with CDS annotation
ExtractSummary.n_with_exons   sequences with EXONS annotation
ExtractSummary.n_no_cds       sequences lacking CDS (skipped for CDS-dependent files)
ExtractSummary.warnings       list of per-sequence warning strings
"""
from __future__ import annotations

import csv
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from lib.fasta_header_parser import HeaderAnnotation, read_fasta_records

log = logging.getLogger('lib.fasta_to_extract')

_NA = 'NA'

# Manifest columns that junctions/architecture/uorf expect
_MANIFEST_HEADER = [
    'transcript_id', 'gene_id', 'region', 'length',
    'source_file', 'seqid', 'strand',
]


# ── Summary ───────────────────────────────────────────────────────────────────

@dataclass
class ExtractSummary:
    n_sequences: int = 0
    n_with_cds:  int = 0
    n_with_exons: int = 0
    n_no_cds:    int = 0
    warnings:    list[str] = field(default_factory=list)


# ── Main entry point ──────────────────────────────────────────────────────────

def build_extract_dir(
    fasta_path: Path,
    extract_dir: Path,
    reference_gff: Optional[Path] = None,
    force: bool = False,
) -> ExtractSummary:
    """
    Read *fasta_path* and write the extract_dir file tree.

    Parameters
    ----------
    fasta_path    : multi-FASTA of full transcripts (UTR + CDS).
    extract_dir   : destination directory (created if absent).
    reference_gff : optional reference GFF3 for WT mode.
                    When supplied, exon/CDS coordinates are taken from the
                    GFF (matched by transcript ID) rather than from FASTA
                    header tags.  Header tags still supply fallback coords
                    if a transcript is absent from the GFF.
    force         : rewrite even if outputs are current.

    Skips writing if all output files are newer than fasta_path unless
    *force* is True.
    """
    extract_dir.mkdir(parents=True, exist_ok=True)

    sentinel = extract_dir / 'manifest.tsv'
    if not force and _is_current(sentinel, fasta_path):
        log.info(f"extract_dir is current ({extract_dir}) — skipping rebuild.")
        summary = ExtractSummary()
        with open(fasta_path) as fh:
            summary.n_sequences = sum(1 for l in fh if l.startswith('>'))
        return summary

    log.info(f"Building extract_dir from FASTA: {fasta_path}"
             + (f" + GFF: {reference_gff}" if reference_gff else " (header-only mode)"))
    records = read_fasta_records(fasta_path)
    summary = ExtractSummary(n_sequences=len(records))

    # Pre-load GFF annotation index when a reference GFF is supplied (WT mode).
    gff_index: dict = {}
    if reference_gff and reference_gff.exists():
        try:
            gff_index = _build_gff_index(reference_gff)
            log.info(f"GFF index built: {len(gff_index)} transcripts")
        except Exception as exc:
            log.warning(f"Could not index GFF ({exc}); falling back to header annotations.")

    manifest_rows: list[dict] = []
    gff_lines: list[str] = []
    utr5_records: list[tuple[str, str]] = []
    cds_records:  list[tuple[str, str]] = []
    utr3_records: list[tuple[str, str]] = []
    mrna_records: list[tuple[str, str]] = []

    gff_lines.append('##gff-version 3')

    for ann, seq in records:
        seq_len = len(seq)
        warns = ann.validate(seq_len)
        for w in warns:
            log.warning(w)
        summary.warnings.extend(warns)

        # ── Merge GFF annotation into the HeaderAnnotation (WT mode) ─────────
        # If a GFF index entry exists for this transcript, its exon/CDS
        # coordinates take precedence over any header-embedded tags.
        if gff_index:
            gff_entry = gff_index.get(ann.seq_id)
            if gff_entry:
                if gff_entry.get('cds'):
                    ann.cds = gff_entry['cds']   # already 0-based half-open
                if gff_entry.get('exons'):
                    ann.exons = gff_entry['exons']
                # Re-derive UTRs now that CDS may have changed
                ann.utr5 = None
                ann.utr3 = None
                ann.derive_utrs(len(seq))
            else:
                log.debug(
                    f"{ann.seq_id}: not found in GFF index; "
                    f"using header annotations."
                )

        if ann.has_cds:
            summary.n_with_cds += 1
        else:
            summary.n_no_cds += 1
            log.warning(
                f"{ann.seq_id}: no CDS annotation (header or GFF) — will appear "
                f"in manifest and mRNA FASTA but CDS/UTR files will be absent."
            )

        if ann.has_exons:
            summary.n_with_exons += 1

        gene_id = ann.gene_id or ann.seq_id

        # ── manifest rows ─────────────────────────────────────────────────────
        # mRNA row (always present)
        manifest_rows.append({
            'transcript_id': ann.seq_id,
            'gene_id':       gene_id,
            'region':        'mRNA',
            'length':        seq_len,
            'source_file':   str(fasta_path),
            'seqid':         ann.seq_id,
            'strand':        ann.strand,
        })
        mrna_records.append((ann.seq_id, seq))

        if ann.has_cds:
            cds_s, cds_e = ann.cds
            cds_seq = seq[cds_s:cds_e]

            # 5' UTR (may be absent for CDS-starting transcripts)
            if ann.utr5 and ann.utr5[1] > ann.utr5[0]:
                u5s, u5e = ann.utr5
                utr5_seq = seq[u5s:u5e]
                manifest_rows.append({
                    'transcript_id': ann.seq_id,
                    'gene_id':       gene_id,
                    'region':        '5UTR',
                    'length':        u5e - u5s,
                    'source_file':   str(fasta_path),
                    'seqid':         ann.seq_id,
                    'strand':        ann.strand,
                })
                utr5_records.append((ann.seq_id, utr5_seq))

            # CDS
            manifest_rows.append({
                'transcript_id': ann.seq_id,
                'gene_id':       gene_id,
                'region':        'CDS',
                'length':        cds_e - cds_s,
                'source_file':   str(fasta_path),
                'seqid':         ann.seq_id,
                'strand':        ann.strand,
            })
            cds_records.append((ann.seq_id, cds_seq))

            # 3' UTR (may be absent)
            if ann.utr3 and ann.utr3[1] > ann.utr3[0]:
                u3s, u3e = ann.utr3
                utr3_seq = seq[u3s:u3e]
                manifest_rows.append({
                    'transcript_id': ann.seq_id,
                    'gene_id':       gene_id,
                    'region':        '3UTR',
                    'length':        u3e - u3s,
                    'source_file':   str(fasta_path),
                    'seqid':         ann.seq_id,
                    'strand':        ann.strand,
                })
                utr3_records.append((ann.seq_id, utr3_seq))

            # ── GFF3 features ─────────────────────────────────────────────────
            gff_lines.extend(
                _gff_for_transcript(ann, seq_len, gene_id)
            )

    # ── Write files ───────────────────────────────────────────────────────────
    _write_manifest(extract_dir / 'manifest.tsv', manifest_rows)
    _write_gff(extract_dir / 'canonical.gff', gff_lines)
    _write_fasta(extract_dir / 'extracted_mRNA.fa',    mrna_records)
    _write_fasta(extract_dir / 'extracted_5UTR.fa',    utr5_records)
    _write_fasta(extract_dir / 'extracted_CDS.fa',     cds_records)
    _write_fasta(extract_dir / 'extracted_3UTR.fa',    utr3_records)

    log.info(
        f"extract_dir built: {summary.n_sequences} sequences, "
        f"{summary.n_with_cds} with CDS, "
        f"{summary.n_with_exons} with EXONS, "
        f"{summary.n_no_cds} lacking CDS annotation."
    )
    return summary


# ── GFF3 synthesis ────────────────────────────────────────────────────────────

def _gff_for_transcript(
    ann: HeaderAnnotation,
    seq_len: int,
    gene_id: str,
) -> list[str]:
    """
    Emit GFF3 lines for one transcript, matching the feature types and
    attribute keys that lib/gff.py and the metric plugins expect:

      gene            (1 per transcript — allows gene_id lookup)
      mRNA            (transcript-level parent)
      exon            (one per exon; derived from EXONS= or inferred from CDS=)
      CDS             (one per CDS segment, tiled with exons if multi-exon)
      five_prime_UTR  (if 5' UTR present)
      three_prime_UTR (if 3' UTR present)

    Coordinates are 1-based inclusive (GFF3 standard).
    seqid = transcript ID (self-contained; gffutils treats each transcript
    as its own "chromosome").
    """
    lines: list[str] = []
    sid   = ann.seq_id
    src   = 'fasta_to_extract'
    strd  = ann.strand
    tx_id = sid

    def gff(feature, start0, end0, attrs):
        # start0/end0 are 0-based half-open → GFF3 1-based inclusive
        return (
            f"{sid}\t{src}\t{feature}\t{start0 + 1}\t{end0}\t"
            f".\t{strd}\t.\t{attrs}"
        )

    # gene
    lines.append(gff('gene', 0, seq_len,
                      f'ID={gene_id};gene_id={gene_id}'))
    # mRNA / transcript
    lines.append(gff('mRNA', 0, seq_len,
                      f'ID={tx_id};Parent={gene_id};'
                      f'transcript_id={tx_id};gene_id={gene_id}'))

    if not ann.has_cds:
        # No CDS: emit a single exon spanning the full transcript
        lines.append(gff('exon', 0, seq_len,
                          f'Parent={tx_id};transcript_id={tx_id}'))
        return lines

    cds_s, cds_e = ann.cds

    # ── Exons ─────────────────────────────────────────────────────────────────
    if ann.has_exons:
        exon_intervals = ann.exons   # already sorted, 0-based half-open
    else:
        # Single-exon transcript spanning the full sequence
        exon_intervals = [(0, seq_len)]

    for es, ee in exon_intervals:
        lines.append(gff('exon', es, ee,
                          f'Parent={tx_id};transcript_id={tx_id}'))

    # ── CDS segments (intersect exons with the CDS range) ────────────────────
    # Each exon that overlaps [cds_s, cds_e) contributes a CDS feature.
    for es, ee in exon_intervals:
        seg_s = max(es, cds_s)
        seg_e = min(ee, cds_e)
        if seg_e > seg_s:
            lines.append(gff('CDS', seg_s, seg_e,
                              f'Parent={tx_id};transcript_id={tx_id}'))

    # ── UTR features ──────────────────────────────────────────────────────────
    if ann.utr5 and ann.utr5[1] > ann.utr5[0]:
        # Split UTR across exons that overlap it
        for es, ee in exon_intervals:
            seg_s = max(es, ann.utr5[0])
            seg_e = min(ee, ann.utr5[1])
            if seg_e > seg_s:
                lines.append(gff('five_prime_UTR', seg_s, seg_e,
                                  f'Parent={tx_id};transcript_id={tx_id}'))

    if ann.utr3 and ann.utr3[1] > ann.utr3[0]:
        for es, ee in exon_intervals:
            seg_s = max(es, ann.utr3[0])
            seg_e = min(ee, ann.utr3[1])
            if seg_e > seg_s:
                lines.append(gff('three_prime_UTR', seg_s, seg_e,
                                  f'Parent={tx_id};transcript_id={tx_id}'))

    # ── start_codon / stop_codon (first/last 3 nt of CDS, genomic coords) ────
    # Needed by junctions.py's reference-point extraction.
    # On + strand: start_codon = cds_s..cds_s+3, stop_codon = cds_e-3..cds_e
    # On − strand: reversed, but since our transcript IS the + strand of the
    # spliced mRNA we always treat it as +.
    lines.append(gff('start_codon', cds_s, cds_s + 3,
                      f'Parent={tx_id};transcript_id={tx_id}'))
    lines.append(gff('stop_codon',  cds_e - 3, cds_e,
                      f'Parent={tx_id};transcript_id={tx_id}'))

    return lines


# ── File writers ──────────────────────────────────────────────────────────────

def _write_manifest(path: Path, rows: list[dict]) -> None:
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_MANIFEST_HEADER, delimiter='\t',
            extrasaction='ignore', lineterminator='\n',
        )
        writer.writeheader()
        writer.writerows(rows)
    log.debug(f"Wrote manifest ({len(rows)} rows) → {path}")


def _write_gff(path: Path, lines: list[str]) -> None:
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    log.debug(f"Wrote canonical.gff ({len(lines)} lines) → {path}")


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, 'w') as fh:
        for seq_id, seq in records:
            fh.write(f'>{seq_id}\n')
            # 60-char wrapped
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + '\n')
    log.debug(f"Wrote {len(records)} records → {path}")


# ── Staleness check ───────────────────────────────────────────────────────────

def _is_current(output: Path, source: Path) -> bool:
    return (
        output.exists()
        and source.exists()
        and output.stat().st_mtime > source.stat().st_mtime
    )


# ── GFF index for WT mode ─────────────────────────────────────────────────────

def _build_gff_index(gff_path: Path) -> dict[str, dict]:
    """
    Parse a reference GFF3 and return a lightweight index:
      {transcript_id: {'cds': (start0, end0), 'exons': [(s0,e0), ...]}}

    Coordinates are 0-based half-open (converted from GFF 1-based inclusive).

    Uses gffutils if available for robust parsing; falls back to a line-by-line
    parser for environments where gffutils is not installed.
    """
    try:
        import gffutils
        return _build_gff_index_gffutils(gff_path)
    except ImportError:
        log.debug("gffutils not available; using line parser for GFF index.")
        return _build_gff_index_lineparser(gff_path)


def _build_gff_index_gffutils(gff_path: Path) -> dict[str, dict]:
    """Build GFF index via gffutils (preferred)."""
    import gffutils, tempfile, os
    db_path = gff_path.with_suffix(gff_path.suffix + '.index.db')
    try:
        db = gffutils.FeatureDB(str(db_path))
    except Exception:
        db = gffutils.create_db(
            str(gff_path), str(db_path),
            merge_strategy='merge', sort_attribute_values=True,
            disable_infer_genes=True, disable_infer_transcripts=True,
            force=True,
        )

    index: dict[str, dict] = {}
    for tx in db.features_of_type(('mRNA', 'transcript')):
        tx_id = tx.id
        exons = sorted(
            db.children(tx, featuretype='exon'),
            key=lambda e: e.start,
        )
        cds_feats = sorted(
            db.children(tx, featuretype='CDS'),
            key=lambda c: c.start,
        )
        entry: dict = {}
        if exons:
            entry['exons'] = [(e.start - 1, e.end) for e in exons]
        if cds_feats:
            entry['cds'] = (
                min(c.start - 1 for c in cds_feats),
                max(c.end       for c in cds_feats),
            )
        if entry:
            index[tx_id] = entry
    return index


def _build_gff_index_lineparser(gff_path: Path) -> dict[str, dict]:
    """Fallback GFF index builder — plain line parser, no dependencies."""
    import re
    import gzip as _gz

    open_fn = _gz.open if str(gff_path).endswith('.gz') else open
    tx_exons: dict[str, list] = {}
    tx_cds:   dict[str, list] = {}

    def _attr(attrs: str, keys) -> Optional[str]:
        for key in keys:
            m = re.search(rf'{key}[= ]"?([^";]+)"?', attrs)
            if m:
                return m.group(1).strip()
        return None

    with open_fn(gff_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            feat   = parts[2]
            start0 = int(parts[3]) - 1
            end0   = int(parts[4])
            attrs  = parts[8]
            tx_id  = _attr(attrs, ('transcript_id', 'Parent', 'Name'))
            if not tx_id:
                continue
            if feat == 'exon':
                tx_exons.setdefault(tx_id, []).append((start0, end0))
            elif feat == 'CDS':
                tx_cds.setdefault(tx_id, []).append((start0, end0))

    index: dict[str, dict] = {}
    all_ids = set(tx_exons) | set(tx_cds)
    for tx_id in all_ids:
        entry: dict = {}
        if tx_id in tx_exons:
            entry['exons'] = sorted(tx_exons[tx_id])
        if tx_id in tx_cds:
            segs = tx_cds[tx_id]
            entry['cds'] = (min(s for s, _ in segs), max(e for _, e in segs))
        index[tx_id] = entry
    return index