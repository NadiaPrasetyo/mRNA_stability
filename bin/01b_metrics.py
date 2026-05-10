#!/usr/bin/env python3
"""
01b_metrics.py
Run lightweight per-transcript metrics defined as plugins under metrics/.

Each metric is a Python module under metrics/<name>.py exposing the plugin
contract documented in metrics/__init__.py.

Outputs land at: $RUNS_ROOT/<dataset>/metrics/<plugin_output>.tsv

Usage:
  ./bin/01b_metrics.py --dataset human_test
  ./bin/01b_metrics.py -d human_test --metric junctions
  ./bin/01b_metrics.py -d human_test --force
  ./bin/01b_metrics.py -d human_test --list-plugins
"""
from __future__ import annotations

import argparse
import importlib
import logging
import os
import sys
from pathlib import Path
from typing import Iterable

# --- Path bootstrap -----------------------------------------------------------
# Import lib/paths.py without requiring the project to be a real package.
_THIS = Path(__file__).resolve()
_PROJECT_ROOT = _THIS.parent.parent
sys.path.insert(0, str(_PROJECT_ROOT))

from lib.paths import PathContext, resolve_paths, add_project_root_to_syspath  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
    datefmt='%H:%M:%S'
)
log = logging.getLogger('metrics')


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
    """True if output exists and is newer than every existing input.

    Missing inputs are skipped (they're the plugin's problem to flag).
    """
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


def _run_one_plugin(name: str, mcfg: dict, paths: PathContext, force: bool) -> bool:
    """Run a single plugin. Returns True on success (including 'skipped'),
    False on failure."""
    plugin_log = logging.getLogger(f'metrics.{name}')
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

    # Staleness check (orchestrator-level; plugin doesn't repeat it)
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


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--dataset', '-d',
                        default=os.environ.get('DATASET'),
                        help="Dataset name (configs/datasets/<name>.yaml)")
    parser.add_argument('--metric', '-m', action='append', default=None,
                        help="Run only this metric. Repeatable. Default: all "
                             "enabled in YAML.")
    parser.add_argument('--force', action='store_true',
                        help="Recompute even if outputs are current.")
    parser.add_argument('--list-plugins', action='store_true',
                        help="List discovered plugin modules and exit.")
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.list_plugins:
        plugins = _list_available_plugins()
        if not plugins:
            print("(no plugins found under metrics/)", file=sys.stderr)
            return
        print("Available plugins:")
        for p in plugins:
            print(f"  {p}")
        return

    if not args.dataset:
        log.error("--dataset (or DATASET env var) required")
        sys.exit(1)

    add_project_root_to_syspath(_PROJECT_ROOT)

    # Resolve paths and load dataset YAML
    try:
        # First pass — no species yet; we read it from YAML and re-resolve.
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

    # Determine which metrics to run
    enabled = _enabled_metrics(config)
    if args.metric:
        # User requested specific metrics; warn if any aren't in enabled set
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

    # Run sequentially. Metrics are cheap; serial is simpler than parallel
    # and lets logs interleave cleanly.
    n_ok = 0
    for name, mcfg in enabled.items():
        if _run_one_plugin(name, mcfg, paths, args.force):
            n_ok += 1

    log.info(f"Metrics complete: {n_ok}/{len(enabled)} succeeded.")
    if n_ok < len(enabled):
        sys.exit(1)


if __name__ == '__main__':
    main()
