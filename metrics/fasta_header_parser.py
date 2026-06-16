"""lib/fasta_header_parser.py
Parse structured annotation tags embedded in FASTA sequence headers.

Supported tags (all 1-based inclusive, matching GFF/GenBank convention):
  CDS=<start>..<end>            e.g.  CDS=101..450
  EXONS=<s1>-<e1>,<s2>-<e2>   e.g.  EXONS=1-200,350-600,800-1000
  UTR5=<start>..<end>           e.g.  UTR5=1..100   (optional)
  UTR3=<start>..<end>           e.g.  UTR3=451..600  (optional)
  GENE=<gene_id>                e.g.  GENE=BRCA1     (optional)
  STRAND=+|-                    e.g.  STRAND=+        (optional; default +)

Tags are whitespace- or pipe-delimited anywhere after the sequence ID:
  >MY_TX|CDS=101..450|EXONS=1-200,350-600
  >MY_TX CDS=101..450 EXONS=1-200,350-600 GENE=BRCA1

All coordinates are stored internally as 0-based half-open Python ranges,
matching gffutils and the rest of the pipeline. Tag coordinates are
converted at parse time so callers never need to think about it.

Public API
──────────
  parse_header(header_line: str) -> HeaderAnnotation
  read_fasta_records(fasta_path: Path) -> list[tuple[HeaderAnnotation, str]]

HeaderAnnotation exposes:
  .seq_id          str
  .gene_id         str  ('' if absent)
  .strand          str  ('+' or '-')
  .cds             tuple[int,int] | None   0-based half-open
  .utr5            tuple[int,int] | None   derived from CDS if absent
  .utr3            tuple[int,int] | None   derived from CDS+len if absent
  .exons           list[tuple[int,int]]    sorted, 0-based half-open
  .has_cds         bool
  .has_exons       bool
  .exon_junctions()             -> list[int]   transcript-coordinate junction positions
  .derive_utrs(seq_len)         fills utr5/utr3 from CDS in place; returns self
  .validate(seq_len)            -> list[str] warning strings
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

Range = Tuple[int, int]   # 0-based half-open


# ── Regexes ───────────────────────────────────────────────────────────────────

_CDS_RE    = re.compile(r'\bCDS=(\d+)\.\.(\d+)\b',    re.IGNORECASE)
_UTR5_RE   = re.compile(r'\bUTR5=(\d+)\.\.(\d+)\b',   re.IGNORECASE)
_UTR3_RE   = re.compile(r'\bUTR3=(\d+)\.\.(\d+)\b',   re.IGNORECASE)
_EXONS_RE  = re.compile(r'\bEXONS=([\d,\-]+)\b',       re.IGNORECASE)
_GENE_RE   = re.compile(r'\bGENE=([^\s|;]+)\b',        re.IGNORECASE)
_STRAND_RE = re.compile(r'\bSTRAND=([+\-])\b',         re.IGNORECASE)
_EXON_INT  = re.compile(r'(\d+)-(\d+)')


def _to0(start1: int, end1: int) -> Range:
    """1-based inclusive → 0-based half-open."""
    return (start1 - 1, end1)


# ── Data class ────────────────────────────────────────────────────────────────

@dataclass
class HeaderAnnotation:
    seq_id:   str
    gene_id:  str = ''
    strand:   str = '+'
    cds:      Optional[Range] = None
    utr5:     Optional[Range] = None
    utr3:     Optional[Range] = None
    exons:    List[Range] = field(default_factory=list)
    raw:      str = ''

    # ── Derived ───────────────────────────────────────────────────────────────

    @property
    def has_cds(self) -> bool:
        return self.cds is not None

    @property
    def has_exons(self) -> bool:
        return bool(self.exons)

    @property
    def cds_start(self) -> Optional[int]:
        return self.cds[0] if self.cds else None

    @property
    def cds_end(self) -> Optional[int]:
        return self.cds[1] if self.cds else None

    # ── Derived UTR fill-in ───────────────────────────────────────────────────

    def derive_utrs(self, seq_len: int) -> 'HeaderAnnotation':
        """Infer UTR ranges from CDS + sequence length if not explicitly set."""
        if self.cds is None:
            return self
        cds_s, cds_e = self.cds
        if self.utr5 is None and cds_s > 0:
            self.utr5 = (0, cds_s)
        if self.utr3 is None and cds_e < seq_len:
            self.utr3 = (cds_e, seq_len)
        return self

    # ── Junction positions in transcript (spliced) coordinates ───────────────

    def exon_junctions(self) -> List[int]:
        """
        Return 0-based transcript positions of each exon-exon junction
        (= cumulative length after each non-final exon).

        This is the position *between* exons in spliced transcript space:
        the first nt of the downstream exon's spliced coordinate.
        Returns [] for 0 or 1 exons.
        """
        if len(self.exons) < 2:
            return []
        junctions: List[int] = []
        cumulative = 0
        for s, e in self.exons[:-1]:
            cumulative += e - s
            junctions.append(cumulative)
        return junctions

    # ── Validation ────────────────────────────────────────────────────────────

    def validate(self, seq_len: int) -> List[str]:
        """Return warning strings for annotation inconsistencies."""
        warns: List[str] = []
        if self.cds:
            s, e = self.cds
            if s < 0 or e > seq_len:
                warns.append(
                    f"{self.seq_id}: CDS {s+1}..{e} out of bounds (seq_len={seq_len})"
                )
            if (e - s) % 3 != 0:
                warns.append(
                    f"{self.seq_id}: CDS length {e - s} nt is not divisible by 3"
                )
        for i, (s, e) in enumerate(self.exons):
            if s < 0 or e > seq_len:
                warns.append(
                    f"{self.seq_id}: exon {i+1} ({s+1}..{e}) out of bounds "
                    f"(seq_len={seq_len})"
                )
        if self.exons:
            total = sum(e - s for s, e in self.exons)
            if total != seq_len:
                warns.append(
                    f"{self.seq_id}: exon total {total} nt ≠ seq len {seq_len} nt. "
                    f"EXONS= must tile the full transcript (UTR + CDS)."
                )
        return warns


# ── Parsing ───────────────────────────────────────────────────────────────────

def parse_header(line: str) -> HeaderAnnotation:
    """Parse a FASTA header line (with or without leading '>') into a
    HeaderAnnotation.  The sequence ID is the first whitespace/pipe token."""
    line = line.lstrip('>')
    seq_id = re.split(r'[\s|]', line, maxsplit=1)[0]

    ann = HeaderAnnotation(seq_id=seq_id, raw=line)

    m = _CDS_RE.search(line)
    if m:
        ann.cds = _to0(int(m.group(1)), int(m.group(2)))

    m = _UTR5_RE.search(line)
    if m:
        ann.utr5 = _to0(int(m.group(1)), int(m.group(2)))

    m = _UTR3_RE.search(line)
    if m:
        ann.utr3 = _to0(int(m.group(1)), int(m.group(2)))

    m = _EXONS_RE.search(line)
    if m:
        ann.exons = sorted(
            [_to0(int(im.group(1)), int(im.group(2)))
             for im in _EXON_INT.finditer(m.group(1))],
            key=lambda r: r[0],
        )

    m = _GENE_RE.search(line)
    if m:
        ann.gene_id = m.group(1)

    m = _STRAND_RE.search(line)
    if m:
        ann.strand = m.group(1)

    return ann


def read_fasta_records(
    fasta_path: Path,
) -> list[tuple[HeaderAnnotation, str]]:
    """Read every record from *fasta_path*.

    Returns a list of (HeaderAnnotation, uppercase_sequence) pairs with
    UTRs derived from CDS + sequence length where not explicitly tagged.
    """
    results: list[tuple[HeaderAnnotation, str]] = []
    current_ann: Optional[HeaderAnnotation] = None
    parts: list[str] = []

    with open(fasta_path, 'r') as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith('>'):
                if current_ann is not None:
                    seq = ''.join(parts).upper().replace('U', 'T')
                    current_ann.derive_utrs(len(seq))
                    results.append((current_ann, seq))
                current_ann = parse_header(line)
                parts = []
            elif line:
                parts.append(line)

    if current_ann is not None:
        seq = ''.join(parts).upper().replace('U', 'T')
        current_ann.derive_utrs(len(seq))
        results.append((current_ann, seq))

    return results