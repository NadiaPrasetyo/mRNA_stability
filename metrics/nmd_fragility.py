"""metrics/nmd_fragility.py
Calculates the "translational fragility" of a transcript by measuring the density 
of out-of-frame stop codons and in-frame 1-nt transition-to-stop codons within 
the NMD-competent zones of the CDS.

Output (<metrics_dir>/nmd_fragility.tsv) columns:
    transcript_id               joinable to manifest.tsv
    gene_id                     joinable across metrics
    strand
    cds_length                  Total length of the extracted CDS (nt)
    nmd_zone_length             Total length of NMD-competent CDS regions (nt)
    n_alt_stops                 Count of out-of-frame TAA/TAG/TGA in the NMD zone
    n_fragile_codons            Count of in-frame CGA/CAG/CAA/TGG in the NMD zone
    alt_stop_density            n_alt_stops / nmd_zone_length
    fragile_codon_density       n_fragile_codons / nmd_zone_length
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterable, Optional

from Bio import SeqIO

from lib.gff import (
    open_or_build_db,
    build_transcript_index,
    lookup_transcript,
    normalise_id,
    split_composite_fasta_id,
)

log = logging.getLogger('metrics.nmd_fragility')

OUTPUT_FILENAME = 'nmd_fragility.tsv'
_NA = 'NA'

_HEADER = [
    'transcript_id', 'gene_id', 'strand', 
    'cds_length', 'nmd_zone_length',
    'n_alt_stops', 'n_fragile_codons',
    'alt_stop_density', 'fragile_codon_density'
]

# Codons that trigger NMD fragility 
ALT_STOPS = {'TAA', 'TAG', 'TGA'}
FRAGILE_CODONS = {'CGA', 'CAG', 'CAA', 'TGG'}


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [
        paths.canonical_gff, 
        paths.manifest_tsv,
        paths.extract_dir / 'extracted_CDS.fa'
    ]


def compute(paths, metric_config, output_path: Path) -> None:
    db = open_or_build_db(paths.canonical_gff,
                          paths.extract_dir / 'canonical.gff.db')

    index = build_transcript_index(db)
    if not index:
        log.error("No transcript features found in canonical.gff. Cannot proceed.")
        _write_tsv(output_path, [])
        return

    fasta_path = paths.extract_dir / 'extracted_CDS.fa'
    if not fasta_path.exists():
        log.error(f"Required FASTA not found: {fasta_path}")
        _write_tsv(output_path, [])
        return

    rows = []
    n_missing, n_skipped = 0, 0

    # We iterate over the FASTA because it provides the physical sequences we need to scan
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Use the updated helper to safely parse the composite header
        parsed = split_composite_fasta_id(record.id, region='CDS')
        if not parsed:
            continue
            
        _, tx_id, _ = parsed

        tx_feature = lookup_transcript(db, tx_id, index)
        if tx_feature is None:
            n_missing += 1
            continue
            
        row = _compute_for_transcript(db, tx_feature, record, manifest_tx_id=tx_id)
        if row is None:
            n_skipped += 1
            continue
            
        rows.append(row)

    if n_missing:
        log.warning(f"{n_missing} transcript(s) in FASTA not found in GFF index.")
    if n_skipped:
        log.warning(f"{n_skipped} transcript(s) skipped (e.g., missing exons/start codons).")

    _write_tsv(output_path, rows)
    log.info(f"Wrote {len(rows)} rows to {output_path.name}")


# ---------------------------------------------------------------------------
# Spliced-coordinate machinery (Ported from junctions.py to keep plugin self-contained)
# ---------------------------------------------------------------------------

def _build_spliced_index(exons, strand: str):
    if strand == '-':
        sorted_exons = sorted(exons, key=lambda e: e.start, reverse=True)
    else:
        sorted_exons = sorted(exons, key=lambda e: e.start)

    spliced_index = []
    cum = 0
    for e in sorted_exons:
        length = e.end - e.start + 1
        spliced_index.append((e.start, e.end, cum, length))
        cum += length

    junctions = []
    cum = 0
    for _, _, _, length in spliced_index[:-1]:
        cum += length
        junctions.append(cum)

    return spliced_index, junctions


def _genomic_to_spliced(genomic_pos: int, spliced_index, strand: str) -> Optional[int]:
    for g_start, g_end, spliced_start, _ in spliced_index:
        if g_start <= genomic_pos <= g_end:
            if strand == '-':
                offset = g_end - genomic_pos
            else:
                offset = genomic_pos - g_start
            return spliced_start + offset
    return None


def _start_codon_genomic_pos(start_codons, cds_features, strand: str) -> Optional[int]:
    if start_codons:
        if strand == '-':
            return max(c.end for c in start_codons)
        return min(c.start for c in start_codons)
    if cds_features:
        if strand == '-':
            return max(c.end for c in cds_features)
        return min(c.start for c in cds_features)
    return None


# ---------------------------------------------------------------------------
# Per-transcript computation
# ---------------------------------------------------------------------------

def _compute_for_transcript(db, transcript, seq_record, manifest_tx_id: str) -> Optional[dict]:
    strand = transcript.strand
    exons = list(db.children(transcript, featuretype='exon'))
    if not exons:
        return None

    spliced_index, junctions = _build_spliced_index(exons, strand)

    cds_features = list(db.children(transcript, featuretype='CDS'))
    start_codons = list(db.children(transcript, featuretype='start_codon'))
    start_g = _start_codon_genomic_pos(start_codons, cds_features, strand)

    start_s = (_genomic_to_spliced(start_g, spliced_index, strand)
               if start_g is not None else None)

    if start_s is None:
        return None

    # Normalise RNA sequences to DNA for standardized codon matching
    cds_seq = str(seq_record.seq).upper().replace('U', 'T')
    cds_len = len(cds_seq)

    nmd_zone_length = 0
    n_alt_stops = 0
    n_fragile_codons = 0

    # 1. Map the NMD-Competent Zones
    # is_competent[i] will be True if the nucleotide at cds_seq[i] 
    # is >50nt upstream of the *next* exon-exon junction.
    is_competent = [False] * cds_len
    
    for i in range(cds_len):
        p = start_s + i
        # Find the first junction downstream of current position 'p'
        next_j = None
        for j in junctions:
            if j > p:
                next_j = j
                break
        
        # The 50nt rule: NMD triggers if distance to junction is strictly > 50
        if next_j is not None and (next_j - p) > 50:
            is_competent[i] = True
            nmd_zone_length += 1

    # 2. Scan for tripwires within the competent zones
    for i in range(cds_len - 2):
        if is_competent[i]:
            codon = cds_seq[i:i+3]
            
            if i % 3 != 0:
                # Out-of-frame: Check for ambush stops
                if codon in ALT_STOPS:
                    n_alt_stops += 1
            else:
                # In-frame: Check for single-transition nonsense mutations
                if codon in FRAGILE_CODONS:
                    n_fragile_codons += 1

    gene_id = ''
    for attr in ('gene_id', 'Parent'):
        if attr in transcript.attributes and transcript.attributes[attr]:
            gene_id = normalise_id(transcript.attributes[attr][0])
            break

    # 3. Calculate Densities
    alt_stop_density = (n_alt_stops / nmd_zone_length) if nmd_zone_length > 0 else None
    fragile_codon_density = (n_fragile_codons / nmd_zone_length) if nmd_zone_length > 0 else None

    return {
        'transcript_id': manifest_tx_id,
        'gene_id': gene_id,
        'strand': strand,
        'cds_length': cds_len,
        'nmd_zone_length': nmd_zone_length,
        'n_alt_stops': n_alt_stops,
        'n_fragile_codons': n_fragile_codons,
        'alt_stop_density': alt_stop_density,
        'fragile_codon_density': fragile_codon_density,
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
