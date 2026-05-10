"""lib/paths.py
Shared path resolver for Python pipeline scripts.

Counterpart to lib/paths.sh: the shell scripts source paths.sh and call
resolve_paths; Python scripts import this module and call resolve_paths.

Both produce the same view of the world for a given dataset.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class PathContext:
    """Resolved paths and dataset metadata for a single pipeline invocation.

    Frozen so plugins can't accidentally mutate it.
    """
    project_root: Path
    runs_root: Path
    dataset: str
    dataset_yaml: Path

    # Dataset-level outputs (shared across tools and metrics)
    run_dir: Path
    extract_dir: Path
    canonical_gff: Path           # filtered GFF, one transcript per gene
    manifest_tsv: Path
    metrics_dir: Path

    # Optional: species name and references root (set if dataset YAML
    # declares 'species:'; consumed by metric plugins that need
    # species-specific reference data e.g. CSC tables).
    species: Optional[str] = None
    references_root: Optional[Path] = None


def project_root_from(script_path: str) -> Path:
    """Project root = parent of the directory containing `script_path`.

    Used by bin/*.py scripts which all live one level under PROJECT_ROOT.
    """
    return Path(script_path).resolve().parent.parent


def resolve_paths(
    dataset: str,
    project_root: Path,
    species: Optional[str] = None,
) -> PathContext:
    """Resolve all paths for a (dataset, project_root) pair.

    Does not source any shell config — this is a pure-Python view used by
    01_extract.py, 01b_metrics.py and metric plugins. The shell side
    (lib/paths.sh) handles cluster + tool configs which Python doesn't need.

    Raises FileNotFoundError if the dataset YAML is missing.
    """
    project_root = Path(project_root).resolve()
    runs_root = Path(os.environ.get('RUNS_ROOT', project_root / 'runs'))

    dataset_yaml = project_root / 'configs' / 'datasets' / f'{dataset}.yaml'
    if not dataset_yaml.is_file():
        raise FileNotFoundError(f"Dataset config not found: {dataset_yaml}")

    run_dir = runs_root / dataset
    extract_dir = run_dir / 'extracted_regions'
    metrics_dir = run_dir / 'metrics'

    references_root: Optional[Path] = None
    if species:
        references_root = project_root / 'data' / 'references' / species

    return PathContext(
        project_root=project_root,
        runs_root=runs_root,
        dataset=dataset,
        dataset_yaml=dataset_yaml,
        run_dir=run_dir,
        extract_dir=extract_dir,
        canonical_gff=extract_dir / 'canonical.gff',
        manifest_tsv=extract_dir / 'manifest.tsv',
        metrics_dir=metrics_dir,
        species=species,
        references_root=references_root,
    )


def add_project_root_to_syspath(project_root: Path) -> None:
    """Allow 'from metrics.<name> import ...' style imports from bin/ scripts."""
    p = str(project_root)
    if p not in sys.path:
        sys.path.insert(0, p)
