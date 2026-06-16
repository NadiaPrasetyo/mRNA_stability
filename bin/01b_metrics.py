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

FASTA mode  (new):
  ./bin/01b_metrics.py --fasta transcripts.fa --species "Homo sapiens"
  ./bin/01b_metrics.py --fasta transcripts.fa --species "Mus musculus" \\
                       --output-dir results/ --metric cai --metric junctions

  • Sequences must be full transcripts (UTR + CDS).
  • Species is used to auto-fetch a GFF3 from NCBI Datasets (primary) or
    Ensembl REST (fallback) for metrics that need annotation
    (junctions, architecture, nmd_fragility_*, uorf).
  • Fetched GFF is cached under <output_dir>/gff_cache/ so subsequent runs
    are instant.

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
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable, NamedTuple, Optional

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
# Plugins NOT in this set are assumed to work from sequence alone.
_GFF_DEPENDENT_PLUGINS: frozenset[str] = frozenset({
    'architecture',
    'junctions',
    'nmd_fragility_core',
    'nmd_fragility_full',
    'nmd_fragility_window',
    'uorf',
})


# ══════════════════════════════════════════════════════════════════════════════
# GFF auto-fetch helpers
# ══════════════════════════════════════════════════════════════════════════════

def _ncbi_taxon_id(species: str) -> Optional[str]:
    """Return NCBI taxonomy ID for *species* (binomial name) via E-utilities."""
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
    """
    Try to download a reference GFF3 from NCBI Datasets for *species*.

    Uses the NCBI Datasets v2 API:
      GET /genome/taxon/{taxon}/annotation_report  → pick latest RefSeq assembly
      GET /genome/accession/{accession}/download?include=GFF3

    The GFF3 is cached at cache_dir/<safe_species>.gff3.gz and the
    decompressed version at cache_dir/<safe_species>.gff3.
    """
    safe = species.replace(' ', '_').lower()
    gff_path = cache_dir / f"{safe}.gff3"
    gz_path = cache_dir / f"{safe}.gff3.gz"

    if gff_path.exists():
        log.info(f"GFF cache hit (NCBI): {gff_path}")
        return gff_path

    taxon_id = _ncbi_taxon_id(species)
    if not taxon_id:
        log.debug(f"No NCBI taxon ID for '{species}'")
        return None

    # Find the latest RefSeq reference assembly accession
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
        log.debug(f"No RefSeq assemblies for taxon {taxon_id}")
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
        # NCBI returns a zip; we need to extract the GFF from it
        import io, zipfile
        req = urllib.request.Request(download_url, headers={'Accept': 'application/zip'})
        with urllib.request.urlopen(req, timeout=120) as r:
            zip_bytes = r.read()
        with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
            gff_members = [n for n in zf.namelist()
                           if n.endswith('.gff') or n.endswith('.gff3')]
            if not gff_members:
                log.debug("No GFF file inside NCBI zip")
                return None
            with zf.open(gff_members[0]) as src, open(gff_path, 'wb') as dst:
                dst.write(src.read())
        log.info(f"NCBI GFF saved → {gff_path}")
        return gff_path
    except Exception as e:
        log.debug(f"NCBI GFF download failed: {e}")
        return None


def _ensembl_species_name(species: str) -> str:
    """Convert 'Homo sapiens' → 'homo_sapiens' for Ensembl REST."""
    return species.strip().lower().replace(' ', '_')


def _fetch_gff_ensembl(species: str, cache_dir: Path) -> Optional[Path]:
    """
    Fallback: download GFF3 from Ensembl REST for *species*.

    Uses:
      GET https://rest.ensembl.org/info/assembly/<species>  → get assembly name
      ftp/https bulk download for GFF3 is too slow; instead we use the
      overlap endpoint to warn users this path is annotation-light, or
      redirect to the Ensembl FTP GFF3.

    Because the Ensembl FTP path is stable and versioned, we build the URL
    from the species name and the current Ensembl release obtained via REST.
    """
    safe = species.replace(' ', '_').lower()
    gff_path = cache_dir / f"{safe}.gff3"

    if gff_path.exists():
        log.info(f"GFF cache hit (Ensembl): {gff_path}")
        return gff_path

    ens_name = _ensembl_species_name(species)

    # Get current Ensembl release number
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

    # Construct FTP URL — Ensembl names the file <Species>.<Assembly>.<release>.gff3.gz
    # We need the assembly name first.
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

    # Capitalise species component correctly for the FTP path
    # e.g. homo_sapiens  →  Homo_sapiens
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
        # Decompress
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
    """
    Return path to a GFF3 for *species*, fetching and caching if needed.

    Strategy:
      1. NCBI Datasets (RefSeq reference assembly — most complete)
      2. Ensembl FTP  (fallback for non-vertebrates / species not in RefSeq)
    """
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
# Synthetic PathContext for FASTA mode
# ══════════════════════════════════════════════════════════════════════════════

def _build_fasta_path_context(
    fasta: Path,
    species: str,
    output_dir: Path,
    gff: Optional[Path],
    engineered: bool = False,
) -> object:
    """
    Construct a PathContext-like object exposing every attribute the real
    metric plugins read.

    Attribute mapping (plugin expectation -> what we provide):
      paths.canonical_gff  -> <extract_dir>/canonical.gff
                              (synthesised from headers, or symlinked WT GFF)
      paths.manifest_tsv   -> <extract_dir>/manifest.tsv
      paths.extract_dir    -> <output_dir>/extract/
      paths.metrics_dir    -> <output_dir>/metrics/
      paths.fasta_path     -> <extract_dir>/extracted_mRNA.fa  (consistent
                              copy so sequence-only plugins share one root)
      paths.species        -> as supplied
      paths.engineered     -> bool flag

    The extract_dir tree is populated by build_extract_dir() before plugins
    run; this object merely holds the resolved paths.
    """
    extract_dir = output_dir / 'extract'

    class _FastaPathContext:
        """Minimal PathContext shim for FASTA-mode runs."""

        def __init__(self):
            # Paths GFF-dependent plugins read (junctions, architecture,
            # uorf, nmd_fragility_*)
            self.canonical_gff: Path = extract_dir / 'canonical.gff'
            self.manifest_tsv:  Path = extract_dir / 'manifest.tsv'
            self.extract_dir:   Path = extract_dir

            # Sequence-only plugins (sequence_basic, cai, codon_aa_counts,
            # stopfree) point here so all plugins share one input root.
            self.fasta_path:    Path = extract_dir / 'extracted_mRNA.fa'

            # Orchestrator / misc
            self.metrics_dir:   Path = output_dir / 'metrics'
            self.output_dir:    Path = output_dir
            self.species:       str  = species
            self.engineered:    bool = engineered
            self.gff_cache_dir: Path = output_dir / 'gff_cache'
            self.dataset_yaml:  Optional[Path] = None
            self.gff_path:      Optional[Path] = gff   # fetched WT GFF or None

            # Mirror any extra fields from the real PathContext so attribute
            # lookups never raise AttributeError.
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
    """Dynamically import metrics.<name>. Raises ImportError if missing."""
    return importlib.import_module(f'metrics.{name}')


def _is_output_current(output_path: Path, input_paths: Iterable[Path]) -> bool:
    """True if output exists and is newer than every existing input."""
    if not output_path.exists():
        return False
    out_mtime = output_path.stat().st_mtime
    for ip in input_paths:
        ip = Path(ip)
        if ip.exists() and ip.stat().st_mtime > out_mtime:
            return False
    return True


def _enabled_metrics(config: dict) -> dict:
    """Return {name: metric_config} for metrics flagged enabled: true."""
    metrics_block = config.get('metrics') or {}
    enabled = {}
    for name, mcfg in metrics_block.items():
        mcfg = mcfg or {}
        if mcfg.get('enabled', False):
            enabled[name] = mcfg
    return enabled


def _list_available_plugins() -> list[str]:
    """Discover plugin modules by listing metrics/*.py."""
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
    paths,           # PathContext or _FastaPathContext
    force: bool,
    gff_available: bool = True,
) -> bool:
    """
    Run a single plugin. Returns True on success (including 'skipped'),
    False on failure.

    In FASTA mode, GFF-dependent plugins are skipped when no GFF is available.
    """
    plugin_log = logging.getLogger(f'metrics.{name}')

    # ── FASTA-mode GFF guard ──────────────────────────────────────────────────
    if not gff_available and name in _GFF_DEPENDENT_PLUGINS:
        plugin_log.warning(
            f"Skipping '{name}': requires GFF annotation which could not be "
            f"fetched for the specified species."
        )
        return True   # not a hard failure — caller asked us to skip gracefully

    try:
        plugin = _load_plugin(name)
    except ImportError as e:
        plugin_log.error(f"Plugin not found: {e}")
        return False

    # Validate plugin contract
    required = ('OUTPUT_FILENAME', 'get_input_paths', 'compute')
    missing = [r for r in required if not hasattr(plugin, r)]
    if missing:
        plugin_log.error(f"Plugin missing attributes: {missing}")
        return False

    output_path = paths.metrics_dir / plugin.OUTPUT_FILENAME
    inputs = list(plugin.get_input_paths(paths, mcfg))

    # Staleness check
    if not force and _is_output_current(output_path, inputs):
        plugin_log.info(f"Output current ({output_path.name}) — skipping.")
        return True

    # Verify inputs exist before launching
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
# FASTA-mode orchestration
# ══════════════════════════════════════════════════════════════════════════════

def _all_plugins_enabled() -> dict:
    """Return {name: {}} for every discovered plugin (FASTA mode default)."""
    return {name: {} for name in _list_available_plugins()
            if not name.startswith('_')}


def _run_fasta_mode(args: argparse.Namespace) -> None:
    """Entry point for --fasta / --species runs."""
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

    # ── Determine which plugins to run ───────────────────────────────────────
    if args.metric:
        requested = {m: {} for m in args.metric}
    else:
        requested = _all_plugins_enabled()

    if not requested:
        log.warning("No plugins to run.")
        return

    # ── GFF handling and extract_dir preparation ─────────────────────────────
    # Both modes end up populating an extract_dir tree that the real plugins
    # can consume unchanged.  The difference is the source of annotation:
    #   engineered: FASTA header tags  (CDS=, EXONS=, …)
    #   WT:         fetched reference GFF + FASTA slicing
    gff_path: Optional[Path] = None
    gff_available: bool = False
    needs_gff_plugins = any(m in _GFF_DEPENDENT_PLUGINS for m in requested)

    if args.engineered:
        log.info(
            "Engineered mode: GFF fetch skipped — extract_dir will be built "
            "from FASTA header annotations (CDS=, EXONS=, UTR5=, UTR3=)."
        )
        gff_available = True   # extract_dir will be valid; no GFF needed
    else:
        if needs_gff_plugins:
            cache_dir = output_dir / 'gff_cache'
            gff_path = resolve_gff(species, cache_dir)
            gff_available = gff_path is not None
            if not gff_available:
                log.warning(
                    "GFF could not be fetched; GFF-dependent metrics will be skipped."
                )

    # ── Build synthetic PathContext ───────────────────────────────────────────
    paths = _build_fasta_path_context(
        fasta, species, output_dir, gff_path,
        engineered=args.engineered,
    )

    # ── Populate extract_dir ─────────────────────────────────────────────────
    # Always run: sequence-only plugins also read extracted_mRNA.fa from here.
    # For WT mode, build_extract_dir uses the fetched GFF if present; for
    # engineered mode it uses header annotations exclusively.
    try:
        from lib.fasta_to_extract import build_extract_dir
        build_extract_dir(
            fasta_path=fasta,
            extract_dir=paths.extract_dir,
            reference_gff=gff_path,        # None in engineered mode
            force=args.force,
        )
    except Exception as exc:
        log.error(f"Failed to build extract_dir: {exc}")
        sys.exit(1)

    # ── Run plugins ──────────────────────────────────────────────────────────
    log.info(f"Running {len(requested)} metric(s): {', '.join(requested)}")
    paths.metrics_dir.mkdir(parents=True, exist_ok=True)

    n_ok = 0
    for name, mcfg in requested.items():
        if _run_one_plugin(name, mcfg, paths, args.force, gff_available):
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
        help="FASTA mode: path to a multi-FASTA of full transcripts (UTR+CDS). "
             "Requires --species.",
    )

    # ── FASTA-mode options ────────────────────────────────────────────────────
    parser.add_argument(
        '--species', '-s',
        metavar='BINOMIAL',
        help="Binomial species name for GFF auto-fetch "
             "(e.g. 'Homo sapiens'). Required in FASTA mode.",
    )
    parser.add_argument(
        '--output-dir', '-o',
        metavar='DIR',
        help="Output directory for FASTA mode. "
             "Defaults to <fasta_dir>/metrics_<fasta_stem>/",
    )
    parser.add_argument(
        '--engineered',
        action='store_true',
        default=False,
        help=(
            "Mark all sequences in the FASTA as engineered/synthetic. "
            "GFF fetch is skipped entirely. Annotation-dependent metrics "
            "(junctions, architecture, nmd_fragility_*, uorf) read their "
            "feature coordinates from FASTA header tags instead:\n"
            "  CDS=<start>..<end>          (1-based inclusive, required for "
            "most GFF-dependent metrics)\n"
            "  EXONS=<s1>-<e1>,<s2>-<e2>  (1-based inclusive, required for "
            "junction/NMD metrics)\n"
            "  UTR5=<start>..<end>         (optional; derived from CDS if absent)\n"
            "  UTR3=<start>..<end>         (optional; derived from CDS if absent)\n"
            "Example header:\n"
            "  >MY_TX|CDS=101..450|EXONS=1-200,350-600"
        ),
    )

    # ── Shared options ────────────────────────────────────────────────────────
    parser.add_argument(
        '--metric', '-m',
        action='append',
        default=None,
        metavar='NAME',
        help="Run only this metric. Repeatable. "
             "Default: all enabled (dataset mode) / all discovered (FASTA mode).",
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

    # ── --list-plugins: works in either mode ──────────────────────────────────
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
    if args.fasta:
        if not args.species:
            parser.error("--species is required when using --fasta")
        if args.engineered and not args.fasta:
            parser.error("--engineered only makes sense with --fasta")
        _run_fasta_mode(args)
    elif args.dataset:
        _run_dataset_mode(args)
    else:
        parser.error("Provide either --fasta FILE --species NAME  or  --dataset NAME")


if __name__ == '__main__':
    main()