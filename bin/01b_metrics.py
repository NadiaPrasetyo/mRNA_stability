#!/usr/bin/env python3
"""
01b_metrics.py
Run lightweight per-transcript metrics defined as plugins under metrics/.

Each metric is a Python module under metrics/<name>.py exposing the plugin
contract documented in metrics/__init__.py.

─────────────────────────────────────────────────────────────────────────────
MODES
─────────────────────────────────────────────────────────────────────────────
Dataset mode  (original behaviour, unchanged):
  ./bin/01b_metrics.py --dataset human_test

FASTA mode  (single-file, legacy):
  ./bin/01b_metrics.py --fasta transcripts.fa --species "Homo sapiens"

Split-FASTA mode  (new):
  ./bin/01b_metrics.py \\
      --fasta-5utr  sequences_5utr.fa \\
      --fasta-cds   baseline_cds.fa   \\
      --fasta-3utr  sequences_3utr.fa \\
      --species "Homo sapiens"

  Rules:
  • --fasta-cds  must contain exactly ONE sequence.  That single CDS is used
    as the shared CDS for every sample.  Junction metrics are therefore
    computed once from this baseline and are identical across all samples.
  • --fasta-5utr / --fasta-3utr must have matching accession IDs (one entry
    per sample).  A missing partner in either file raises an error unless
    --allow-missing-utrs is set, in which case the absent region is treated
    as zero-length.
  • Full per-sample transcript is assembled in memory as:
        5UTR_seq + CDS_seq + 3UTR_seq
    with FASTA header tags synthesised automatically so downstream plugins
    see the same annotated layout they would from a single-FASTA run.
  • --engineered works identically to single-FASTA mode (no GFF fetch).
  • --metric and --force work the same as in all other modes.

Shared flags:
  ./bin/01b_metrics.py ... --metric junctions   # run only this plugin
  ./bin/01b_metrics.py ... --force              # ignore staleness check
  ./bin/01b_metrics.py ... --list-plugins       # enumerate metrics/ modules
─────────────────────────────────────────────────────────────────────────────
"""
from __future__ import annotations

import argparse
import importlib
import logging
import os
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

# ── Path bootstrap ────────────────────────────────────────────────────────────
_THIS = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from lib.paths import PathContext, resolve_paths, add_project_root_to_syspath  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('metrics')

# ── Metrics that require GFF annotation ───────────────────────────────────────
_GFF_DEPENDENT_PLUGINS: frozenset[str] = frozenset({
    'architecture',
    'junctions',
    'nmd_fragility_core',
    'nmd_fragility_full',
    'nmd_fragility_window',
    'uorf',
})


# ══════════════════════════════════════════════════════════════════════════════
# GFF auto-fetch helpers  (unchanged)
# ══════════════════════════════════════════════════════════════════════════════

def _ncbi_taxon_id(species: str) -> Optional[str]:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=taxonomy&term={urllib.parse.quote(species)}&retmode=json"
    )
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            import json
            data = json.loads(r.read())
        ids = data.get('esearchresult', {}).get('idlist', [])
        return ids[0] if ids else None
    except Exception as e:
        log.debug(f"NCBI taxonomy lookup failed: {e}")
        return None


def _fetch_gff_ncbi(species: str, cache_dir: Path) -> Optional[Path]:
    safe = species.replace(' ', '_').lower()
    gff_path = cache_dir / f"{safe}.gff3"

    if gff_path.exists():
        log.info(f"GFF cache hit (NCBI): {gff_path}")
        return gff_path

    taxon_id = _ncbi_taxon_id(species)
    if not taxon_id:
        log.debug(f"No NCBI taxon ID for '{species}'")
        return None

    report_url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon_id}"
        f"/annotation_report?page_size=5&filters.reference_only=true"
        f"&filters.assembly_source=refseq"
    )
    try:
        req = urllib.request.Request(report_url, headers={'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=20) as r:
            import json
            report = json.loads(r.read())
    except Exception as e:
        log.debug(f"NCBI annotation report failed: {e}")
        return None

    reports = report.get('reports', [])
    if not reports:
        return None

    accession = reports[0].get('accession')
    if not accession:
        return None

    log.info(f"Downloading GFF3 from NCBI for {species} (assembly {accession}) …")
    download_url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}"
        f"/download?include_annotation_type=GENOME_GFF"
    )
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        import io, zipfile
        req = urllib.request.Request(download_url, headers={'Accept': 'application/zip'})
        with urllib.request.urlopen(req, timeout=120) as r:
            zip_bytes = r.read()
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
            gff_members = [n for n in zf.namelist()
                           if n.endswith('.gff') or n.endswith('.gff3')]
            if not gff_members:
                return None
            with zf.open(gff_members[0]) as src, open(gff_path, 'wb') as dst:
                dst.write(src.read())
        log.info(f"NCBI GFF saved → {gff_path}")
        return gff_path
    except Exception as e:
        log.debug(f"NCBI GFF download failed: {e}")
        return None


def _ensembl_species_name(species: str) -> str:
    return species.strip().lower().replace(' ', '_')


def _fetch_gff_ensembl(species: str, cache_dir: Path) -> Optional[Path]:
    safe = species.replace(' ', '_').lower()
    gff_path = cache_dir / f"{safe}.gff3"

    if gff_path.exists():
        log.info(f"GFF cache hit (Ensembl): {gff_path}")
        return gff_path

    ens_name = _ensembl_species_name(species)
    try:
        req = urllib.request.Request(
            'https://rest.ensembl.org/info/software',
            headers={'Content-Type': 'application/json'},
        )
        with urllib.request.urlopen(req, timeout=15) as r:
            import json
            info = json.loads(r.read())
        release = info.get('release', 'current')
    except Exception as e:
        log.debug(f"Ensembl release lookup failed: {e}")
        release = 'current'

    try:
        req = urllib.request.Request(
            f'https://rest.ensembl.org/info/assembly/{ens_name}?content-type=application/json',
            headers={'Content-Type': 'application/json'},
        )
        with urllib.request.urlopen(req, timeout=15) as r:
            assembly_info = json.loads(r.read())
        assembly = assembly_info.get('assembly_name', '')
    except Exception as e:
        log.debug(f"Ensembl assembly lookup failed: {e}")
        return None

    if not assembly:
        return None

    ftp_species = ens_name[0].upper() + ens_name[1:]
    gff_url = (
        f"https://ftp.ensembl.org/pub/release-{release}/gff3/{ens_name}/"
        f"{ftp_species}.{assembly}.{release}.gff3.gz"
    )
    gz_path = cache_dir / f"{safe}.gff3.gz"

    log.info(f"Downloading GFF3 from Ensembl for {species} (release {release}) …")
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(gff_url, gz_path)
        import gzip, shutil
        with gzip.open(gz_path, 'rb') as f_in, open(gff_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_path.unlink(missing_ok=True)
        log.info(f"Ensembl GFF saved → {gff_path}")
        return gff_path
    except Exception as e:
        log.debug(f"Ensembl GFF download failed: {e}")
        gz_path.unlink(missing_ok=True)
        return None


def resolve_gff(species: str, cache_dir: Path) -> Optional[Path]:
    gff = _fetch_gff_ncbi(species, cache_dir)
    if gff:
        return gff
    log.info("NCBI fetch unsuccessful — trying Ensembl …")
    gff = _fetch_gff_ensembl(species, cache_dir)
    if gff:
        return gff
    log.warning(
        f"Could not obtain a GFF3 for '{species}'. "
        "GFF-dependent metrics will be skipped."
    )
    return None


# ══════════════════════════════════════════════════════════════════════════════
# Split-FASTA assembly helpers
# ══════════════════════════════════════════════════════════════════════════════

def _read_raw_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file and return {seq_id: uppercase_sequence}."""
    seqs: dict[str, str] = {}
    current_id: Optional[str] = None
    parts: list[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(parts).upper().replace('U', 'T')
                current_id = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)

    if current_id is not None:
        seqs[current_id] = ''.join(parts).upper().replace('U', 'T')

    return seqs


def _read_single_fasta_sequence(path: Path) -> Tuple[str, str]:
    """
    Read a FASTA expected to contain exactly one record.
    Returns (seq_id, sequence).  Raises ValueError if 0 or >1 records found.
    """
    seqs = _read_raw_fasta(path)
    if len(seqs) == 0:
        raise ValueError(f"CDS FASTA is empty: {path}")
    if len(seqs) > 1:
        raise ValueError(
            f"CDS FASTA must contain exactly one sequence, "
            f"found {len(seqs)}: {path}\n"
            f"  IDs: {', '.join(list(seqs)[:5])}{'…' if len(seqs) > 5 else ''}"
        )
    return next(iter(seqs.items()))


def assemble_split_fastas(
    path_5utr: Path,
    path_cds: Path,
    path_3utr: Path,
    allow_missing_utrs: bool = False,
) -> List[Tuple[str, str, int, int, int]]:
    """
    Combine per-sample 5UTR / shared CDS / per-sample 3UTR into full
    transcript records.

    Returns a list of tuples:
      (sample_id, full_sequence, utr5_len, cds_len, utr3_len)

    The caller uses utr5_len and cds_len to synthesise the FASTA header tags
    (CDS=, UTR5=, UTR3=) that fasta_to_extract.py expects.

    Parameters
    ----------
    path_5utr            : FASTA of 5' UTR sequences, one per sample.
    path_cds             : FASTA with a SINGLE CDS sequence (shared baseline).
    path_3utr            : FASTA of 3' UTR sequences, one per sample.
                           Accession IDs must match path_5utr.
    allow_missing_utrs   : if True, samples missing a UTR partner get a
                           zero-length UTR rather than raising an error.
    """
    cds_id, cds_seq = _read_single_fasta_sequence(path_cds)
    cds_len = len(cds_seq)
    log.info(f"Shared CDS: '{cds_id}'  length={cds_len} nt")

    seqs_5 = _read_raw_fasta(path_5utr)
    seqs_3 = _read_raw_fasta(path_3utr)

    ids_5 = set(seqs_5)
    ids_3 = set(seqs_3)

    if ids_5 != ids_3:
        only_5 = ids_5 - ids_3
        only_3 = ids_3 - ids_5
        msg_parts = []
        if only_5:
            msg_parts.append(f"only in 5UTR: {sorted(only_5)[:5]}")
        if only_3:
            msg_parts.append(f"only in 3UTR: {sorted(only_3)[:5]}")
        mismatch_msg = "; ".join(msg_parts)
        if not allow_missing_utrs:
            raise ValueError(
                f"5UTR and 3UTR FASTAs have mismatched accessions. "
                f"{mismatch_msg}\n"
                f"Use --allow-missing-utrs to treat absent regions as "
                f"zero-length."
            )
        log.warning(
            f"5UTR / 3UTR accession mismatch ({mismatch_msg}). "
            "Missing UTRs will be zero-length."
        )

    all_ids = sorted(ids_5 | ids_3)
    records: List[Tuple[str, str, int, int, int]] = []

    for sample_id in all_ids:
        utr5_seq = seqs_5.get(sample_id, '')
        utr3_seq = seqs_3.get(sample_id, '')
        full_seq = utr5_seq + cds_seq + utr3_seq

        utr5_len = len(utr5_seq)
        utr3_len = len(utr3_seq)

        records.append((sample_id, full_seq, utr5_len, cds_len, utr3_len))

    log.info(
        f"Assembled {len(records)} sample transcript(s) from split FASTAs "
        f"({len(seqs_5)} 5UTR, 1 CDS, {len(seqs_3)} 3UTR)."
    )
    return records


def _write_assembled_fasta(
    records: List[Tuple[str, str, int, int, int]],
    out_path: Path,
) -> None:
    """
    Write assembled transcripts to *out_path* with embedded annotation tags.

    Header format:
      ><sample_id>|CDS=<cds_start>..<cds_end>|UTR5=<u5s>..<u5e>|UTR3=<u3s>..<u3e>

    All coordinates are 1-based inclusive (matching GFF / the parser convention).
    UTR5= / UTR3= tags are omitted when the UTR is zero-length.
    """
    with open(out_path, 'w') as fh:
        for sample_id, seq, utr5_len, cds_len, utr3_len in records:
            # 1-based inclusive coordinates
            cds_start = utr5_len + 1
            cds_end   = utr5_len + cds_len

            tags = [f"CDS={cds_start}..{cds_end}"]
            if utr5_len > 0:
                tags.append(f"UTR5=1..{utr5_len}")
            if utr3_len > 0:
                u3_start = cds_end + 1
                u3_end   = cds_end + utr3_len
                tags.append(f"UTR3={u3_start}..{u3_end}")

            header = '|'.join([sample_id] + tags)
            fh.write(f'>{header}\n')
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + '\n')

    log.info(f"Assembled FASTA written → {out_path}  ({len(records)} records)")


# ══════════════════════════════════════════════════════════════════════════════
# Shared CDS junction FASTA helper
# ══════════════════════════════════════════════════════════════════════════════

def _write_cds_only_fasta(
    cds_path: Path,
    out_path: Path,
) -> None:
    """
    Write a single-record annotated FASTA containing only the shared CDS.

    The header tags set CDS to span the entire sequence (no UTRs), which
    means junction metrics see only the CDS exon structure.  When the CDS
    FASTA has no EXONS= tag in its header the full CDS is treated as a
    single exon (which is the common case for a baseline CDS input).

    This file is written to extract_dir/cds_junctions.fa and consumed by
    a dedicated junctions run if --junctions-from-cds is active (default
    in split-FASTA mode).
    """
    cds_id, cds_seq = _read_single_fasta_sequence(cds_path)
    cds_len = len(cds_seq)
    header  = f"{cds_id}|CDS=1..{cds_len}"

    with open(out_path, 'w') as fh:
        fh.write(f'>{header}\n')
        for i in range(0, len(cds_seq), 60):
            fh.write(cds_seq[i:i+60] + '\n')

    log.info(f"CDS-only junction FASTA written → {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# Synthetic PathContext (shared)
# ══════════════════════════════════════════════════════════════════════════════

def _build_fasta_path_context(
    fasta: Path,
    species: str,
    output_dir: Path,
    gff: Optional[Path],
    engineered: bool = False,
) -> object:
    extract_dir = output_dir / 'extract'

    class _FastaPathContext:
        def __init__(self):
            self.canonical_gff: Path = extract_dir / 'canonical.gff'
            self.manifest_tsv:  Path = extract_dir / 'manifest.tsv'
            self.extract_dir:   Path = extract_dir
            self.fasta_path:    Path = extract_dir / 'extracted_mRNA.fa'
            self.metrics_dir:   Path = output_dir / 'metrics'
            self.output_dir:    Path = output_dir
            self.species:       str  = species
            self.engineered:    bool = engineered
            self.gff_cache_dir: Path = output_dir / 'gff_cache'
            self.dataset_yaml:  Optional[Path] = None
            self.gff_path:      Optional[Path] = gff

            try:
                import dataclasses as _dc
                _flds = [f.name for f in _dc.fields(PathContext)]
            except (TypeError, Exception):
                _flds = list(getattr(PathContext, '_fields', []))
            for _f in _flds:
                if not hasattr(self, _f):
                    setattr(self, _f, None)

    return _FastaPathContext()


# ══════════════════════════════════════════════════════════════════════════════
# Original helpers (unchanged)
# ══════════════════════════════════════════════════════════════════════════════

def _load_yaml(path: Path) -> dict:
    try:
        import yaml
    except ImportError:
        log.error("PyYAML not installed (pip install pyyaml).")
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f) or {}


def _load_plugin(name: str):
    return importlib.import_module(f'metrics.{name}')


def _is_output_current(output_path: Path, input_paths: Iterable[Path]) -> bool:
    if not output_path.exists():
        return False
    out_mtime = output_path.stat().st_mtime
    for ip in input_paths:
        ip = Path(ip)
        if ip.exists() and ip.stat().st_mtime > out_mtime:
            return False
    return True


def _enabled_metrics(config: dict) -> dict:
    metrics_block = config.get('metrics') or {}
    enabled = {}
    for name, mcfg in metrics_block.items():
        mcfg = mcfg or {}
        if mcfg.get('enabled', False):
            enabled[name] = mcfg
    return enabled


def _list_available_plugins() -> list[str]:
    metrics_dir = _PROJECT_ROOT / 'metrics'
    if not metrics_dir.is_dir():
        return []
    return sorted(
        p.stem for p in metrics_dir.glob('*.py')
        if p.stem != '__init__' and not p.stem.startswith('_')
    )


def _run_one_plugin(
    name: str,
    mcfg: dict,
    paths,
    force: bool,
    gff_available: bool = True,
) -> bool:
    plugin_log = logging.getLogger(f'metrics.{name}')

    if not gff_available and name in _GFF_DEPENDENT_PLUGINS:
        plugin_log.warning(
            f"Skipping '{name}': requires GFF annotation which could not be "
            f"fetched for the specified species."
        )
        return True

    try:
        plugin = _load_plugin(name)
    except ImportError as e:
        plugin_log.error(f"Plugin not found: {e}")
        return False

    required = ('OUTPUT_FILENAME', 'get_input_paths', 'compute')
    missing = [r for r in required if not hasattr(plugin, r)]
    if missing:
        plugin_log.error(f"Plugin missing attributes: {missing}")
        return False

    output_path = paths.metrics_dir / plugin.OUTPUT_FILENAME
    inputs = list(plugin.get_input_paths(paths, mcfg))

    if not force and _is_output_current(output_path, inputs):
        plugin_log.info(f"Output current ({output_path.name}) — skipping.")
        return True

    missing_inputs = [str(p) for p in inputs if not Path(p).exists()]
    if missing_inputs:
        plugin_log.error(
            "Required input(s) missing:\n  " + "\n  ".join(missing_inputs))
        return False

    paths.metrics_dir.mkdir(parents=True, exist_ok=True)
    plugin_log.info(f"Computing → {output_path}")
    try:
        plugin.compute(paths, mcfg, output_path)
    except Exception as e:
        plugin_log.exception(f"Plugin raised: {e}")
        return False

    if not output_path.exists():
        plugin_log.error(f"Plugin completed but did not produce {output_path}")
        return False

    plugin_log.info(f"Done ({output_path.stat().st_size} bytes)")
    return True


# ══════════════════════════════════════════════════════════════════════════════
# FASTA-mode orchestration  (single-file, unchanged behaviour)
# ══════════════════════════════════════════════════════════════════════════════

def _all_plugins_enabled() -> dict:
    return {name: {} for name in _list_available_plugins()
            if not name.startswith('_')}


def _run_fasta_mode(args: argparse.Namespace) -> None:
    """Entry point for --fasta / --species runs (single combined FASTA)."""
    fasta = Path(args.fasta).expanduser().resolve()
    if not fasta.exists():
        log.error(f"FASTA not found: {fasta}")
        sys.exit(1)

    species: str = args.species
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir \
        else fasta.parent / f"metrics_{fasta.stem}"
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"FASTA mode — file      : {fasta}")
    log.info(f"            species    : {species}")
    log.info(f"            output     : {output_dir}")
    log.info(f"            engineered : {args.engineered}")

    _run_fasta_core(
        fasta=fasta,
        species=species,
        output_dir=output_dir,
        engineered=args.engineered,
        metric=args.metric,
        force=args.force,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Split-FASTA mode orchestration  (NEW)
# ══════════════════════════════════════════════════════════════════════════════

def _run_split_fasta_mode(args: argparse.Namespace) -> None:
    """
    Entry point for --fasta-5utr / --fasta-cds / --fasta-3utr runs.

    Steps
    ─────
    1. Validate the three input paths.
    2. Read and validate the CDS (must be exactly one sequence).
    3. Read 5UTR and 3UTR FASTAs; check accession parity.
    4. Assemble per-sample full transcripts (5UTR + CDS + 3UTR) with
       embedded CDS= / UTR5= / UTR3= header tags.
    5. Write the assembled FASTA to <output_dir>/assembled_transcripts.fa.
    6. Delegate to the shared FASTA core (which runs GFF fetch + plugins).

    Junctions are calculated from the shared CDS, so every sample has
    identical junction metrics — this is intentional and reflects the
    design constraint that only UTR sequences vary across samples.
    """
    path_5utr = Path(args.fasta_5utr).expanduser().resolve()
    path_cds  = Path(args.fasta_cds).expanduser().resolve()
    path_3utr = Path(args.fasta_3utr).expanduser().resolve()

    for label, p in [('5UTR', path_5utr), ('CDS', path_cds), ('3UTR', path_3utr)]:
        if not p.exists():
            log.error(f"{label} FASTA not found: {p}")
            sys.exit(1)

    species: str = args.species
    output_dir = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else path_cds.parent / f"metrics_{path_cds.stem}_split"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Split-FASTA mode")
    log.info(f"  5UTR FASTA  : {path_5utr}")
    log.info(f"  CDS  FASTA  : {path_cds}  (shared baseline — must be 1 sequence)")
    log.info(f"  3UTR FASTA  : {path_3utr}")
    log.info(f"  Species     : {species}")
    log.info(f"  Output      : {output_dir}")
    log.info(f"  Engineered  : {args.engineered}")
    log.info(
        "  Note: junction metrics are computed from the shared CDS only "
        "and will be identical across all samples."
    )

    # ── Assemble full transcripts ─────────────────────────────────────────────
    try:
        assembled = assemble_split_fastas(
            path_5utr=path_5utr,
            path_cds=path_cds,
            path_3utr=path_3utr,
            allow_missing_utrs=getattr(args, 'allow_missing_utrs', False),
        )
    except ValueError as exc:
        log.error(str(exc))
        sys.exit(1)

    assembled_fasta = output_dir / 'assembled_transcripts.fa'
    _write_assembled_fasta(assembled, assembled_fasta)

    # ── Write the CDS-only FASTA for junction reference ───────────────────────
    # Plugins (junctions, nmd_fragility_*) that read from extract_dir will use
    # the assembled transcripts, but we log explicitly that the CDS-derived
    # junction structure is shared.  No separate run is needed; the CDS= tag
    # in each assembled header ensures all samples map to the same CDS coords
    # (offset by their individual 5UTR length).  Users can inspect the shared
    # CDS structure via the extract_dir/extracted_CDS.fa output.
    cds_ref_fasta = output_dir / 'cds_reference.fa'
    _write_cds_only_fasta(path_cds, cds_ref_fasta)
    log.info(
        f"CDS reference (junction baseline) written → {cds_ref_fasta}. "
        "The CDS structure is identical for all samples."
    )

    # ── Delegate to shared FASTA core ────────────────────────────────────────
    _run_fasta_core(
        fasta=assembled_fasta,
        species=species,
        output_dir=output_dir,
        engineered=args.engineered,
        metric=args.metric,
        force=args.force,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Shared FASTA core  (used by both single-FASTA and split-FASTA modes)
# ══════════════════════════════════════════════════════════════════════════════

def _run_fasta_core(
    fasta: Path,
    species: str,
    output_dir: Path,
    engineered: bool,
    metric: Optional[list],
    force: bool,
) -> None:
    """
    Given a fully assembled + annotated FASTA, run GFF fetch (if needed)
    and execute all requested metric plugins.

    This is the shared implementation consumed by both _run_fasta_mode and
    _run_split_fasta_mode so the two entry points stay DRY.
    """
    if metric:
        requested = {m: {} for m in metric}
    else:
        requested = _all_plugins_enabled()

    if not requested:
        log.warning("No plugins to run.")
        return

    # ── GFF handling ─────────────────────────────────────────────────────────
    gff_path: Optional[Path] = None
    gff_available: bool = False
    needs_gff_plugins = any(m in _GFF_DEPENDENT_PLUGINS for m in requested)

    if engineered:
        log.info(
            "Engineered mode: GFF fetch skipped — extract_dir will be built "
            "from FASTA header annotations (CDS=, EXONS=, UTR5=, UTR3=)."
        )
        gff_available = True
    else:
        if needs_gff_plugins:
            cache_dir = output_dir / 'gff_cache'
            gff_path = resolve_gff(species, cache_dir)
            gff_available = gff_path is not None
            if not gff_available:
                log.warning(
                    "GFF could not be fetched; GFF-dependent metrics will be skipped."
                )

    # ── Build PathContext ─────────────────────────────────────────────────────
    paths = _build_fasta_path_context(
        fasta, species, output_dir, gff_path,
        engineered=engineered,
    )

    # ── Populate extract_dir ─────────────────────────────────────────────────
    try:
        from lib.fasta_to_extract import build_extract_dir
        build_extract_dir(
            fasta_path=fasta,
            extract_dir=paths.extract_dir,
            reference_gff=gff_path,
            force=force,
        )
    except Exception as exc:
        log.error(f"Failed to build extract_dir: {exc}")
        sys.exit(1)

    # ── Run plugins ──────────────────────────────────────────────────────────
    log.info(f"Running {len(requested)} metric(s): {', '.join(requested)}")
    paths.metrics_dir.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    for name, mcfg in requested.items():
        if _run_one_plugin(name, mcfg, paths, force, gff_available):
            n_ok += 1

    log.info(f"Metrics complete: {n_ok}/{len(requested)} succeeded.")
    if n_ok < len(requested):
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# Dataset-mode orchestration (original, unchanged)
# ══════════════════════════════════════════════════════════════════════════════

def _run_dataset_mode(args: argparse.Namespace) -> None:
    add_project_root_to_syspath(_PROJECT_ROOT)

    try:
        paths = resolve_paths(args.dataset, _PROJECT_ROOT)
    except FileNotFoundError as e:
        log.error(str(e))
        sys.exit(1)

    config = _load_yaml(paths.dataset_yaml)
    species = config.get('species')
    if species:
        paths = resolve_paths(args.dataset, _PROJECT_ROOT, species=species)

    log.info(f"Dataset: {args.dataset}"
             + (f" (species: {species})" if species else ""))

    enabled = _enabled_metrics(config)
    if args.metric:
        requested = {}
        for m in args.metric:
            if m in enabled:
                requested[m] = enabled[m]
            elif m in (config.get('metrics') or {}):
                log.warning(f"Metric '{m}' is in YAML but not enabled — running anyway.")
                requested[m] = (config['metrics'][m] or {})
            else:
                log.warning(f"Metric '{m}' has no config block — running with empty config.")
                requested[m] = {}
        enabled = requested

    if not enabled:
        log.warning("No metrics enabled. Add a 'metrics:' block to your dataset YAML.")
        return

    log.info(f"Running {len(enabled)} metric(s): {', '.join(enabled)}")
    paths.metrics_dir.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    for name, mcfg in enabled.items():
        if _run_one_plugin(name, mcfg, paths, args.force):
            n_ok += 1

    log.info(f"Metrics complete: {n_ok}/{len(enabled)} succeeded.")
    if n_ok < len(enabled):
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Mode selection ────────────────────────────────────────────────────────
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        '--dataset', '-d',
        default=os.environ.get('DATASET'),
        metavar='NAME',
        help="Dataset mode: name of configs/datasets/<name>.yaml",
    )
    mode.add_argument(
        '--fasta', '-f',
        metavar='FILE',
        help=(
            "Single-FASTA mode: path to a multi-FASTA of full transcripts "
            "(UTR+CDS) with embedded header tags.  Requires --species.  "
            "For separate UTR/CDS inputs use --fasta-5utr/--fasta-cds/--fasta-3utr."
        ),
    )
    mode.add_argument(
        '--fasta-5utr',
        metavar='FILE',
        dest='fasta_5utr',
        help=(
            "Split-FASTA mode: FASTA of 5' UTR sequences, one per sample.  "
            "Must be combined with --fasta-cds and --fasta-3utr."
        ),
    )

    # ── Split-FASTA companions (not in the mutex group — validated manually) ──
    parser.add_argument(
        '--fasta-cds',
        metavar='FILE',
        dest='fasta_cds',
        help=(
            "Split-FASTA mode: FASTA containing exactly ONE CDS sequence "
            "(the shared baseline).  Junction metrics are derived from this "
            "sequence and are identical for all samples."
        ),
    )
    parser.add_argument(
        '--fasta-3utr',
        metavar='FILE',
        dest='fasta_3utr',
        help=(
            "Split-FASTA mode: FASTA of 3' UTR sequences, one per sample.  "
            "Accession IDs must match --fasta-5utr."
        ),
    )
    parser.add_argument(
        '--allow-missing-utrs',
        action='store_true',
        default=False,
        dest='allow_missing_utrs',
        help=(
            "Split-FASTA mode: if a sample has a 5UTR entry but no matching "
            "3UTR (or vice versa), treat the absent region as zero-length "
            "instead of raising an error."
        ),
    )

    # ── Common FASTA / split-FASTA options ───────────────────────────────────
    parser.add_argument(
        '--species', '-s',
        metavar='BINOMIAL',
        help="Binomial species name for GFF auto-fetch "
             "(e.g. 'Homo sapiens'). Required in FASTA and split-FASTA modes.",
    )
    parser.add_argument(
        '--output-dir', '-o',
        metavar='DIR',
        help=(
            "Output directory.  "
            "Defaults to <fasta_dir>/metrics_<fasta_stem>/ (single-FASTA) or "
            "<cds_dir>/metrics_<cds_stem>_split/ (split-FASTA)."
        ),
    )
    parser.add_argument(
        '--engineered',
        action='store_true',
        default=False,
        help=(
            "Mark sequences as engineered/synthetic.  GFF fetch is skipped.  "
            "In split-FASTA mode, annotation tags (CDS=, UTR5=, UTR3=) are "
            "synthesised automatically from the UTR/CDS lengths; --engineered "
            "is therefore implied and this flag is accepted for explicitness.\n"
            "In single-FASTA mode, header tags must be supplied manually:\n"
            "  CDS=<start>..<end>          (1-based inclusive)\n"
            "  EXONS=<s1>-<e1>,<s2>-<e2>  (1-based inclusive)\n"
            "  UTR5=<start>..<end>         (optional)\n"
            "  UTR3=<start>..<end>         (optional)"
        ),
    )

    # ── Shared options ────────────────────────────────────────────────────────
    parser.add_argument(
        '--metric', '-m',
        action='append',
        default=None,
        metavar='NAME',
        help="Run only this metric. Repeatable. "
             "Default: all enabled (dataset mode) / all discovered (FASTA modes).",
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help="Recompute even if outputs are current.",
    )
    parser.add_argument(
        '--list-plugins',
        action='store_true',
        help="List discovered plugin modules and exit.",
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # ── --list-plugins ────────────────────────────────────────────────────────
    if args.list_plugins:
        plugins = _list_available_plugins()
        if not plugins:
            print("(no plugins found under metrics/)", file=sys.stderr)
            return
        print("Available plugins:")
        print(f"  {'Plugin':<30}  {'Needs GFF (WT)':<18}  {'Engineered header tags'}")
        print(f"  {'-'*30}  {'-'*18}  {'-'*40}")
        _tag_info = {
            'architecture':         'CDS=',
            'junctions':            'CDS=, EXONS=',
            'nmd_fragility_core':   'CDS=, EXONS=',
            'nmd_fragility_full':   'CDS=, EXONS=',
            'nmd_fragility_window': 'CDS=, EXONS=',
            'uorf':                 'CDS=',
        }
        for p in plugins:
            needs_gff = p in _GFF_DEPENDENT_PLUGINS
            tags = _tag_info.get(p, '—')
            gff_marker = 'yes' if needs_gff else 'no'
            print(f"  {p:<30}  {gff_marker:<18}  {tags}")
        return

    # ── Route to the correct mode ─────────────────────────────────────────────
    if args.fasta_5utr:
        # Split-FASTA mode: all three parts required
        missing_args = []
        if not args.fasta_cds:
            missing_args.append('--fasta-cds')
        if not args.fasta_3utr:
            missing_args.append('--fasta-3utr')
        if not args.species:
            missing_args.append('--species')
        if missing_args:
            parser.error(
                f"Split-FASTA mode requires: {', '.join(missing_args)}"
            )
        # In split-FASTA mode, header tags are auto-generated, so engineered
        # is always effectively true for the assembled FASTA.  Force the flag
        # so fasta_to_extract uses header annotations rather than GFF lookup.
        args.engineered = True
        _run_split_fasta_mode(args)

    elif args.fasta:
        if not args.species:
            parser.error("--species is required when using --fasta")
        _run_fasta_mode(args)

    elif args.dataset:
        _run_dataset_mode(args)

    else:
        parser.error(
            "Provide one of:\n"
            "  --fasta FILE --species NAME              (single-FASTA mode)\n"
            "  --fasta-5utr F --fasta-cds F --fasta-3utr F --species NAME  "
            "(split-FASTA mode)\n"
            "  --dataset NAME                           (dataset mode)"
        )


if __name__ == '__main__':
    main()