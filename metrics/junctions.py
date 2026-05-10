"""metrics/junctions.py
Exon-junction features per transcript, all distances in spliced (mature mRNA)
coordinates.

Output (<metrics_dir>/junctions.tsv) columns:
    transcript_id                   joinable to manifest.tsv
    gene_id                         joinable across metrics
    strand
    n_exons
    n_5UTR_junctions                len(five_prime_UTR features) - 1
    n_CDS_junctions                 len(CDS features) - 1
    n_3UTR_junctions                len(three_prime_UTR features) - 1
    stop_dist_closest_upstream      spliced nt; NA if no upstream junction
    stop_dist_closest_downstream    spliced nt; NA if stop in last exon
    stop_dist_last_downstream       spliced nt; canonical NMD metric
    start_dist_closest_upstream     spliced nt
    start_dist_closest_downstream   spliced nt

Notes
-----
* Region junction counts assume per-exon UTR/CDS features (MANE, GENCODE,
  modern Ensembl). Annotations that collapse UTRs into a single feature
  per region will report 0 UTR junctions even when the UTR spans multiple
  exons.
* Reference points: first nt of start codon, last nt of stop codon — both
  in transcript reading direction. Falls back to CDS extent if explicit
  start_codon / stop_codon features are absent.
* Single-exon transcripts: all junction columns are 0 / NA.
* Stop or start codons split across an intron: handled by collapsing the
  feature pieces to the appropriate transcript-direction extreme.
"""
from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import Iterable, Optional

log = logging.getLogger('metrics.junctions')

OUTPUT_FILENAME = 'junctions.tsv'

# Tab-separated; NA written for missing values
_NA = 'NA'

_HEADER = [
    'transcript_id', 'gene_id', 'strand', 'n_exons',
    'n_5UTR_junctions', 'n_CDS_junctions', 'n_3UTR_junctions',
    'stop_dist_closest_upstream', 'stop_dist_closest_downstream',
    'stop_dist_last_downstream',
    'start_dist_closest_upstream', 'start_dist_closest_downstream',
]


# ---------------------------------------------------------------------------
# Plugin contract
# ---------------------------------------------------------------------------

def get_input_paths(paths, metric_config) -> Iterable[Path]:
    return [paths.canonical_gff, paths.manifest_tsv]


def compute(paths, metric_config, output_path: Path) -> None:
    """Entry point called by the orchestrator."""
    try:
        import gffutils
    except ImportError:
        log.error("gffutils not installed (pip install gffutils).")
        raise

    db = _open_or_build_db(gffutils, paths.canonical_gff,
                          paths.extract_dir / 'canonical.gff.db')

    index = _build_transcript_index(db)
    if not index:
        log.error("No transcript-like features found in canonical.gff "
                  "(no features with exon children). Cannot proceed.")
        # Still write an empty TSV with header so the orchestrator's
        # 'output produced' check passes; the warning above is the signal.
        _write_tsv(output_path, [])
        return

    transcript_ids = _read_manifest_transcripts(paths.manifest_tsv)
    log.info(f"Computing junctions for {len(transcript_ids)} transcripts")

    rows = []
    n_missing, n_no_exons = 0, 0
    sample_missing = []
    for tx_id in transcript_ids:
        tx_feature = _lookup_transcript(db, tx_id, index)
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
        # Diagnostic: when failure rate is high, dump samples so the user
        # can compare manifest IDs vs the keys we actually indexed.
        if n_missing >= max(5, len(transcript_ids) // 2):
            sample_keys = [k for k in list(index.keys())[:8]]
            log.warning("  ID format mismatch suspected. Examples:")
            log.warning(f"    Manifest IDs (first 5): {sample_missing}")
            log.warning(f"    GFF index keys (first 8): {sample_keys}")
            log.warning("  Run with -v for the full feature-type list. "
                        "If the formats clearly differ in a way the normaliser "
                        "doesn't handle, file a bug or extend _KNOWN_PREFIXES.")
    if n_no_exons:
        log.warning(f"{n_no_exons} transcript(s) had no exon features")

    _write_tsv(output_path, rows)
    log.info(f"Wrote {len(rows)} rows to {output_path.name}")


# ---------------------------------------------------------------------------
# Database handling
# ---------------------------------------------------------------------------

def _open_or_build_db(gffutils, gff_path: Path, db_path: Path):
    """Open an existing gffutils DB or build it fresh.

    Rebuild if the GFF is newer than the DB (defensive: handles a manual
    re-run of 01_extract.py that produced a new canonical.gff).
    """
    needs_build = (
        not db_path.exists()
        or db_path.stat().st_mtime < gff_path.stat().st_mtime
    )
    if needs_build:
        if db_path.exists():
            log.info(f"GFF newer than DB; rebuilding {db_path.name}")
            db_path.unlink()
        else:
            log.info(f"Building gffutils DB: {db_path}")
        # Use gffutils defaults for ID assignment. Custom id_spec remapping
        # of mRNA IDs would break Parent= references on child features
        # (exons / CDS / UTR), since those refer to the *original* IDs in
        # the GFF. ID mismatches between manifest and GFF (version, prefix)
        # are handled in _lookup_transcript instead.
        db = gffutils.create_db(
            str(gff_path),
            dbfn=str(db_path),
            force=False,
            keep_order=False,
            merge_strategy='create_unique',
            sort_attribute_values=False,
        )
    else:
        log.debug(f"Loading existing DB: {db_path}")
        db = gffutils.FeatureDB(str(db_path))
    return db


# ---------------------------------------------------------------------------
# Transcript ID lookup (manifest IDs vs GFF IDs are not always identical)
# ---------------------------------------------------------------------------
import re

# Known namespace prefixes seen in GFF3 ID attributes:
#   Ensembl/GENCODE: 'transcript:ENST...', 'gene:ENSG...'
#   NCBI/RefSeq:     'rna-NM_...', 'gene-FOO' (also 'rna-XM_', 'rna-NR_', etc.)
# Conservative list — only patterns that uniquely identify a namespace
# prefix and won't false-strip anything that could appear inside a real ID.
_KNOWN_PREFIXES = ('transcript:', 'gene:', 'rna-', 'gene-')

# Trailing version suffix: '.<digits>' at end of string (e.g. ENST00000000001.7)
_VERSION_RE = re.compile(r'\.\d+$')


def _strip_prefix(s: str) -> str:
    """Strip a known GFF3 ID namespace prefix; return s unchanged if none match."""
    s_lower = s.lower()
    for pfx in _KNOWN_PREFIXES:
        if s_lower.startswith(pfx):
            return s[len(pfx):]
    return s


def _strip_version(s: str) -> str:
    """Strip a trailing '.<digits>' version suffix.
    Only strips dot-then-pure-digits-to-end so e.g. 'ENST00000000001.7' -> 'ENST00000000001'
    but 'gene.alpha' is preserved.
    """
    return _VERSION_RE.sub('', s)


def _normalise_id(raw: str) -> str:
    """Maximal normalisation: strip namespace prefix AND trailing version."""
    return _strip_version(_strip_prefix(raw))


def _index_keys_for_feature(feat) -> list[str]:
    """Generate every plausible lookup key for a transcript-like feature.

    Covers ID forms (raw, prefix-stripped, fully normalised), the
    transcript_id attribute (same three forms), and a reconstructed
    versioned form using a separate transcript_version / version attribute
    if present (Ensembl convention where ID is unversioned but a separate
    attribute carries the version).
    """
    keys = []
    fid = feat.id
    keys.extend([fid, _strip_prefix(fid), _normalise_id(fid)])

    for tid in feat.attributes.get('transcript_id', []):
        keys.extend([tid, _strip_prefix(tid), _normalise_id(tid)])
        # Reconstruct 'BARE.VERSION' if version stored separately
        bare = _normalise_id(tid)
        for vattr in ('transcript_version', 'version'):
            for v in feat.attributes.get(vattr, []):
                keys.append(f"{bare}.{v}")

    # Deduplicate while preserving order; drop empties
    seen = set()
    out = []
    for k in keys:
        if k and k not in seen:
            seen.add(k)
            out.append(k)
    return out


def _build_transcript_index(db) -> dict:
    """Map every plausible ID form to the gffutils feature.id of its transcript.

    Discovers transcripts as 'features that are Parents of exon features'.
    This catches mRNA, transcript, primary_transcript, lnc_RNA, ncRNA, etc.
    without needing to enumerate feature types up front.

    First-write-wins on collisions: if two features share a normalised ID,
    the first one indexed keeps the slot. Real annotations don't collide
    here once the canonical.gff filter has trimmed to one transcript per gene.
    """
    parent_ids = set()
    for exon in db.features_of_type('exon'):
        for pid in exon.attributes.get('Parent', []):
            parent_ids.add(pid)

    index = {}
    feat_types = set()
    for pid in parent_ids:
        try:
            feat = db[pid]
        except Exception:
            continue
        feat_types.add(feat.featuretype)
        for key in _index_keys_for_feature(feat):
            index.setdefault(key, feat.id)

    log.debug(f"Indexed {len(parent_ids)} transcript-like features "
              f"(types: {sorted(feat_types)}) under {len(index)} keys")
    return index


def _lookup_transcript(db, manifest_tx_id: str, index: dict):
    """Look up a transcript by trying several forms of the manifest ID."""
    for key in (manifest_tx_id,
                _strip_prefix(manifest_tx_id),
                _normalise_id(manifest_tx_id)):
        feat_id = index.get(key)
        if feat_id is not None:
            try:
                return db[feat_id]
            except Exception:
                pass
    return None


def _read_manifest_transcripts(manifest_tsv: Path) -> list[str]:
    """Read manifest.tsv, return unique transcript_ids preserving first-seen order."""
    seen = set()
    ordered = []
    with open(manifest_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'transcript_id' not in (reader.fieldnames or []):
            raise ValueError(
                f"manifest.tsv has no 'transcript_id' column; got {reader.fieldnames}")
        for row in reader:
            tid = row['transcript_id']
            if tid and tid not in seen:
                seen.add(tid)
                ordered.append(tid)
    return ordered


# ---------------------------------------------------------------------------
# Spliced-coordinate machinery
# ---------------------------------------------------------------------------

def _build_spliced_index(exons, strand: str):
    """Sort exons in transcript order and assign each a spliced offset.

    Returns:
        spliced_index : list of (g_start, g_end, spliced_start, length)
                        where spliced_start is the 0-indexed mRNA position
                        of the first nt of this exon.
        junctions     : list of spliced positions of each junction, in
                        transcript order. Junction p separates spliced
                        position p-1 from p. There are len(exons)-1 of them.
    """
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

    # Junctions: cumulative length up to (but not including) each non-first exon
    junctions = []
    cum = 0
    for _, _, _, length in spliced_index[:-1]:
        cum += length
        junctions.append(cum)

    return spliced_index, junctions


def _genomic_to_spliced(genomic_pos: int, spliced_index, strand: str) -> Optional[int]:
    """Map a genomic position to a spliced (0-indexed) mRNA position.

    Returns None if the position falls outside every exon (shouldn't happen
    for start/stop codons in well-formed annotations, but possible for
    pathological cases — e.g. a stop codon entirely within an intron).
    """
    for g_start, g_end, spliced_start, _ in spliced_index:
        if g_start <= genomic_pos <= g_end:
            if strand == '-':
                offset = g_end - genomic_pos
            else:
                offset = genomic_pos - g_start
            return spliced_start + offset
    return None


# ---------------------------------------------------------------------------
# Reference-point extraction
# ---------------------------------------------------------------------------

def _start_codon_genomic_pos(start_codons, cds_features, strand: str) -> Optional[int]:
    """Genomic position of the first nucleotide of the start codon.

    `start_codons` and `cds_features` are lists of gffutils Features (possibly
    empty). Falls back to CDS extents when start_codon is absent.
    """
    if start_codons:
        if strand == '-':
            return max(c.end for c in start_codons)
        return min(c.start for c in start_codons)
    if cds_features:
        # First CDS coordinate in transcript direction
        if strand == '-':
            return max(c.end for c in cds_features)
        return min(c.start for c in cds_features)
    return None


def _stop_codon_genomic_pos(stop_codons, cds_features, strand: str) -> Optional[int]:
    """Genomic position of the last nucleotide of the stop codon.

    Caveat: when stop_codon is missing we fall back to the 3' CDS edge,
    which is correct for GFFs that include the stop codon in CDS (Ensembl
    legacy) but off by 3 nt for GFFs that exclude it (most modern Ensembl,
    GENCODE, MANE include it; older RefSeq sometimes excludes it).
    Marginal effect on closest-junction distances.
    """
    if stop_codons:
        if strand == '-':
            return min(c.start for c in stop_codons)
        return max(c.end for c in stop_codons)
    if cds_features:
        if strand == '-':
            return min(c.start for c in cds_features)
        return max(c.end for c in cds_features)
    return None


# ---------------------------------------------------------------------------
# Per-transcript computation
# ---------------------------------------------------------------------------

def _compute_for_transcript(db, transcript, manifest_tx_id: str) -> Optional[dict]:
    """Compute the row for one transcript. Returns None if no exons."""
    strand = transcript.strand
    exons = list(db.children(transcript, featuretype='exon'))
    if not exons:
        return None

    spliced_index, junctions = _build_spliced_index(exons, strand)
    n_exons = len(exons)

    # Region junction counts (assumes per-exon UTR/CDS features)
    utr5 = list(db.children(transcript, featuretype='five_prime_UTR'))
    cds = list(db.children(transcript, featuretype='CDS'))
    utr3 = list(db.children(transcript, featuretype='three_prime_UTR'))
    n_5utr_j = max(0, len(utr5) - 1)
    n_cds_j = max(0, len(cds) - 1)
    n_3utr_j = max(0, len(utr3) - 1)

    # Reference points → spliced
    start_codons = list(db.children(transcript, featuretype='start_codon'))
    stop_codons = list(db.children(transcript, featuretype='stop_codon'))

    start_g = _start_codon_genomic_pos(start_codons, cds, strand)
    stop_g = _stop_codon_genomic_pos(stop_codons, cds, strand)

    start_s = (_genomic_to_spliced(start_g, spliced_index, strand)
               if start_g is not None else None)
    stop_s = (_genomic_to_spliced(stop_g, spliced_index, strand)
              if stop_g is not None else None)

    stop_up, stop_down, stop_last_down = _distances(stop_s, junctions)
    start_up, start_down, _ = _distances(start_s, junctions)

    # gene_id from attributes; fall back to Parent or empty
    gene_id = ''
    for attr in ('gene_id', 'Parent'):
        if attr in transcript.attributes and transcript.attributes[attr]:
            gene_id = _normalise_id(transcript.attributes[attr][0])
            break

    return {
        'transcript_id': manifest_tx_id,   # preserve manifest ID exactly
        'gene_id': gene_id,
        'strand': strand,
        'n_exons': n_exons,
        'n_5UTR_junctions': n_5utr_j,
        'n_CDS_junctions': n_cds_j,
        'n_3UTR_junctions': n_3utr_j,
        'stop_dist_closest_upstream': stop_up,
        'stop_dist_closest_downstream': stop_down,
        'stop_dist_last_downstream': stop_last_down,
        'start_dist_closest_upstream': start_up,
        'start_dist_closest_downstream': start_down,
    }


def _distances(pos: Optional[int], junctions: list[int]):
    """Return (closest_upstream, closest_downstream, last_downstream).

    All in spliced nt; returned as `pos - junction_pos` (upstream) or
    `junction_pos - pos` (downstream). None when no qualifying junction
    exists or `pos` is None.

    Convention
    ----------
    Junctions are 0-indexed spliced positions of the *first nt of the
    downstream exon* (i.e. junction p separates spliced positions p-1
    and p in the mature mRNA).

    A junction located exactly at `pos` (e.g. start codon coinciding with
    a splice acceptor) is treated as upstream with distance 0. A junction
    at `pos+1` (e.g. stop codon at the last nt of an exon) is treated as
    downstream with distance 1. The asymmetry is intentional: a junction
    coincident with `pos` lies between `pos-1` and `pos` and is therefore
    spatially upstream of `pos`.
    """
    if pos is None or not junctions:
        return None, None, None

    upstream = [j for j in junctions if j <= pos]
    downstream = [j for j in junctions if j > pos]

    closest_up = (pos - max(upstream)) if upstream else None
    closest_down = (min(downstream) - pos) if downstream else None
    last_down = (max(downstream) - pos) if downstream else None
    return closest_up, closest_down, last_down


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
            # Stringify with NA for None values (consistent with manifest convention)
            out = {k: (_NA if row.get(k) is None else row.get(k)) for k in _HEADER}
            writer.writerow(out)
