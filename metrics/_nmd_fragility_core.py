"""metrics/_nmd_fragility_core.py
Shared implementation for nmd_fragility plugins.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterable, Optional

from lib.gff import (
    open_or_build_db,
    build_transcript_index,
    lookup_transcript,
    normalise_id,
    split_composite_fasta_id,
)

_NA = 'NA'
_NMD_THRESHOLD = 50      # nt; strict (> 50, not >= 50). Default for core/full
                         # and the fallback default for window's
                         # configurable `nmd_threshold` key.
_VALID_MODELS = ('core', 'full', 'window')

_HEADER = [
    'transcript_id', 'gene_id', 'strand', 'model',
    'cds_length', 'zone_length',
    'n_transition_fragile_codons',
    'n_transversion_fragile_codons',
    'n_snv_fragile_codons',
    'n_alt_stop_codons',
    'transition_fragile_codon_density',
    'transversion_fragile_codon_density',
    'snv_fragile_codon_density',
    'alt_stop_codon_density',
    'transition_fraction_of_snv_fragile',
]

ALT_STOP_CODONS = {'TAA', 'TAG', 'TGA'}

# Sense codons one transition away from a stop codon.
TRANSITION_FRAGILE_CODONS = {'CGA', 'CAG', 'CAA', 'TGG'}

# Sense codons one transversion away from a stop codon.
TRANSVERSION_FRAGILE_CODONS = {
    'AAA', 'AAG',
    'AGA', 'GGA',
    'GAA', 'GAG',
    'TAC', 'TAT',
    'TCA', 'TCG',
    'TGC', 'TGT',
    'TTA', 'TTG',
}

SNV_FRAGILE_CODONS = TRANSITION_FRAGILE_CODONS | TRANSVERSION_FRAGILE_CODONS


# ---------------------------------------------------------------------------
# Entry points used by the variant plugins
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [
        paths.canonical_gff,
        paths.manifest_tsv,
        paths.extract_dir / 'extracted_CDS.fa',
    ]

def run(paths, metric_config, output_path: Path, model: str, log: logging.Logger) -> None:
    if model not in _VALID_MODELS:
        raise ValueError(f"model must be one of {_VALID_MODELS}, got {model!r}")

    db = open_or_build_db(paths.canonical_gff, paths.extract_dir / 'canonical.gff.db')

    index = build_transcript_index(db)
    if not index:
        log.error("No transcript features found in canonical.gff. Cannot proceed.")
        _write_tsv(output_path, [])
        return

    fasta_path = paths.extract_dir / 'extracted_CDS.fa'

    rows: list[dict] = []
    n_unparsed, n_missing, n_skipped = 0, 0, 0

    for record_id, sequence in _iter_fasta(fasta_path):
        parsed = split_composite_fasta_id(record_id, region='CDS')
        if not parsed:
            n_unparsed += 1
            continue
        _, tx_id, _ = parsed

        tx_feature = lookup_transcript(db, tx_id, index)
        if tx_feature is None:
            n_missing += 1
            continue

        row = _compute_for_transcript(
            db, tx_feature, sequence,
            manifest_tx_id=tx_id, model=model, metric_config=metric_config
        )
        if row is None:
            n_skipped += 1
            continue

        rows.append(row)

    if n_unparsed:
        log.warning(f"{n_unparsed} FASTA header(s) did not match the composite pattern.")
    if n_missing:
        log.warning(f"{n_missing} transcript(s) in FASTA not found in GFF index.")
    if n_skipped:
        log.warning(f"{n_skipped} transcript(s) skipped.")

    _write_tsv(output_path, rows)
    log.info(f"Wrote {len(rows)} rows to {output_path.name} (model={model}).")


# ---------------------------------------------------------------------------
# FASTA reading
# ---------------------------------------------------------------------------

def _iter_fasta(path: Path):
    header: Optional[str] = None
    chunks: list[str] = []
    with open(path) as f:
        for raw in f:
            line = raw.rstrip('\r\n')
            if not line: continue
            if line.startswith('>'):
                if header is not None: yield header, ''.join(chunks)
                header = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None: yield header, ''.join(chunks)


# ---------------------------------------------------------------------------
# Spliced-coordinate machinery
# ---------------------------------------------------------------------------

def _build_spliced_index(exons, strand: str):
    if strand == '-': sorted_exons = sorted(exons, key=lambda e: e.start, reverse=True)
    else: sorted_exons = sorted(exons, key=lambda e: e.start)

    spliced_index, junctions, cum = [], [], 0
    for e in sorted_exons:
        length = e.end - e.start + 1
        spliced_index.append((e.start, e.end, cum, length))
        cum += length

    cum = 0
    for _, _, _, length in spliced_index[:-1]:
        cum += length
        junctions.append(cum)

    return spliced_index, junctions

def _genomic_to_spliced(genomic_pos: int, spliced_index, strand: str) -> Optional[int]:
    for g_start, g_end, spliced_start, _ in spliced_index:
        if g_start <= genomic_pos <= g_end:
            if strand == '-': offset = g_end - genomic_pos
            else: offset = genomic_pos - g_start
            return spliced_start + offset
    return None

def _start_codon_genomic_pos(start_codons, cds_features, strand: str) -> Optional[int]:
    if start_codons:
        if strand == '-': return max(c.end for c in start_codons)
        return min(c.start for c in start_codons)
    if cds_features:
        if strand == '-': return max(c.end for c in cds_features)
        return min(c.start for c in cds_features)
    return None


# ---------------------------------------------------------------------------
# Competent-zone definitions
# ---------------------------------------------------------------------------

def _competent_zone_core(start_s: int, cds_len: int, junctions: list[int], metric_config: dict) -> list[bool]:
    is_competent = [False] * cds_len
    ji = 0
    while ji < len(junctions) and junctions[ji] <= start_s: ji += 1
    for i in range(cds_len):
        p = start_s + i
        while ji < len(junctions) and junctions[ji] <= p: ji += 1
        if ji < len(junctions) and (junctions[ji] - p) > _NMD_THRESHOLD:
            is_competent[i] = True
    return is_competent

def _competent_zone_full(start_s: int, cds_len: int, junctions: list[int], metric_config: dict) -> list[bool]:
    is_competent = [False] * cds_len
    if not junctions: return is_competent
    last_junction = junctions[-1]
    zone_end = max(0, min(cds_len, last_junction - start_s - _NMD_THRESHOLD))
    for i in range(zone_end): is_competent[i] = True
    return is_competent

def _competent_zone_window(start_s: int, cds_len: int, junctions: list[int], metric_config: dict) -> list[bool]:
    """Mark positions inside a configurable window anchored to the last junction.

    Config keys (all optional):
        window_size    (int, default 120): length of the proximal scanning window in nt.
        nmd_threshold  (int, default 50):  size of the NMD-immune safe zone in nt;
                                           ignored when apply_nmd_rule is False.
        apply_nmd_rule (bool, default True): if False, the window runs right up to
                                             the junction (exclusive); the position
                                             at the junction is still excluded.

    If the available upstream CDS is shorter than the configured window, the
    zone is clamped to the CDS boundary and `zone_length` reflects the
    actual scanned length (so densities normalise correctly).
    """
    is_competent = [False] * cds_len
    if not junctions:
        return is_competent

    # Pull config with explicit type coercion (guards against YAML int-as-float etc.)
    window_size    = int(metric_config.get('window_size', 120))
    nmd_threshold  = int(metric_config.get('nmd_threshold', _NMD_THRESHOLD))
    apply_nmd_rule = bool(metric_config.get('apply_nmd_rule', True))

    if window_size <= 0:
        raise ValueError(f"window_size must be a positive int, got {window_size}")
    if nmd_threshold < 0:
        raise ValueError(f"nmd_threshold must be non-negative, got {nmd_threshold}")

    target_j = junctions[-1]

    # Window end in spliced coordinates (exclusive upper bound).
    if apply_nmd_rule:
        global_end = target_j - nmd_threshold
    else:
        global_end = target_j

    global_start = global_end - window_size

    # Translate to 0-indexed CDS array coordinates with symmetric clamping.
    cds_start = max(0, min(cds_len, global_start - start_s))
    cds_end   = max(0, min(cds_len, global_end   - start_s))

    for i in range(cds_start, cds_end):
        is_competent[i] = True

    return is_competent


_MODEL_FN = {
    'core': _competent_zone_core,
    'full': _competent_zone_full,
    'window': _competent_zone_window,
}


# ---------------------------------------------------------------------------
# Per-transcript computation
# ---------------------------------------------------------------------------

def _compute_for_transcript(db, transcript, sequence: str, manifest_tx_id: str, model: str, metric_config: dict) -> Optional[dict]:
    strand = transcript.strand
    exons = list(db.children(transcript, featuretype='exon'))
    if not exons: return None

    spliced_index, junctions = _build_spliced_index(exons, strand)

    cds_features = list(db.children(transcript, featuretype='CDS'))
    start_codons = list(db.children(transcript, featuretype='start_codon'))
    start_g = _start_codon_genomic_pos(start_codons, cds_features, strand)
    start_s = (_genomic_to_spliced(start_g, spliced_index, strand) if start_g is not None else None)
    if start_s is None: return None

    cds_seq = sequence.upper().replace('U', 'T')
    cds_len = len(cds_seq)

    # Pass the metric_config to the model function
    is_competent = _MODEL_FN[model](start_s, cds_len, junctions, metric_config)
    zone_length = sum(is_competent)

    n_transition_fragile_codons = 0
    n_transversion_fragile_codons = 0
    n_snv_fragile_codons = 0
    n_alt_stop_codons = 0

    for i in range(cds_len - 2):
        # For NMD purposes, the biologically relevant position is the stop
        # codon's termination point, not merely the first base of the codon.
        # This matches the junction-distance convention used elsewhere in the
        # pipeline: stop-codon distances are measured from the last nt of the
        # stop codon in transcript coordinates.
        if not is_competent[i + 2]:
            continue

        codon = cds_seq[i:i + 3]
        if i % 3 != 0:
            if codon in ALT_STOP_CODONS:
                n_alt_stop_codons += 1
        else:
            is_transition_fragile = codon in TRANSITION_FRAGILE_CODONS
            is_transversion_fragile = codon in TRANSVERSION_FRAGILE_CODONS

            if is_transition_fragile:
                n_transition_fragile_codons += 1
            if is_transversion_fragile:
                n_transversion_fragile_codons += 1
            if is_transition_fragile or is_transversion_fragile:
                n_snv_fragile_codons += 1

    if zone_length > 0:
        transition_fragile_codon_density = n_transition_fragile_codons / zone_length
        transversion_fragile_codon_density = n_transversion_fragile_codons / zone_length
        snv_fragile_codon_density = n_snv_fragile_codons / zone_length
        alt_stop_codon_density = n_alt_stop_codons / zone_length
    else:
        transition_fragile_codon_density = None
        transversion_fragile_codon_density = None
        snv_fragile_codon_density = None
        alt_stop_codon_density = None

    if n_snv_fragile_codons > 0:
        transition_fraction_of_snv_fragile = n_transition_fragile_codons / n_snv_fragile_codons
    else:
        transition_fraction_of_snv_fragile = None

    gene_id = ''
    for attr in ('gene_id', 'Parent'):
        if attr in transcript.attributes and transcript.attributes[attr]:
            gene_id = normalise_id(transcript.attributes[attr][0])
            break

    return {
        'transcript_id': manifest_tx_id, 'gene_id': gene_id, 'strand': strand,
        'model': model, 'cds_length': cds_len, 'zone_length': zone_length,
        'n_transition_fragile_codons': n_transition_fragile_codons,
        'n_transversion_fragile_codons': n_transversion_fragile_codons,
        'n_snv_fragile_codons': n_snv_fragile_codons,
        'n_alt_stop_codons': n_alt_stop_codons,
        'transition_fragile_codon_density': transition_fragile_codon_density,
        'transversion_fragile_codon_density': transversion_fragile_codon_density,
        'snv_fragile_codon_density': snv_fragile_codon_density,
        'alt_stop_codon_density': alt_stop_codon_density,
        'transition_fraction_of_snv_fragile': transition_fraction_of_snv_fragile,
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _write_tsv(output_path: Path, rows: list[dict]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=_HEADER, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        for row in rows:
            out = {k: (_NA if row.get(k) is None else row.get(k)) for k in _HEADER}
            writer.writerow(out)
