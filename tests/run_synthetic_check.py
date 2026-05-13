"""Run the NMD-fragility sanity check outside pytest.

Useful when a test fails and you want to inspect the actual TSV output
directly, or when iterating on the plugin code and want a fast pass /
fail signal without the pytest harness.

Usage (from the repo root):

    python tests/run_synthetic_check.py
    python tests/run_synthetic_check.py --keep-outputs /tmp/nmd_check

Exit code 0 if every assertion holds, 1 otherwise.
"""
from __future__ import annotations

import argparse
import csv
import logging
import shutil
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from metrics._nmd_fragility_core import run as run_core  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures"

# Strand → GFF filename. The FASTA, manifest, and expected outputs are
# identical between strands by design.
GFF_BY_STRAND = {
    "plus":  "synthetic_5x200.gff3",
    "minus": "synthetic_5x200_minus.gff3",
}

# (model, metric_config, expected_dict)
CASES = [
    ("nearest",       {},                            {"nmd_zone_length": 600, "n_fragile_codons": 4, "n_alt_stops": 2}),
    ("any",           {},                            {"nmd_zone_length": 750, "n_fragile_codons": 5, "n_alt_stops": 2}),
    ("distal_window", {},                            {"nmd_zone_length": 120, "n_fragile_codons": 3, "n_alt_stops": 1}),
    ("distal_window", {"apply_nmd_rule": False},     {"nmd_zone_length": 120, "n_fragile_codons": 3, "n_alt_stops": 0}),
]


def stage_fixtures(workdir: Path, strand: str) -> SimpleNamespace:
    extract = workdir / "extract"
    extract.mkdir(parents=True, exist_ok=True)
    shutil.copy(FIXTURE_DIR / GFF_BY_STRAND[strand],     extract / "canonical.gff")
    shutil.copy(FIXTURE_DIR / "synthetic_5x200_CDS.fa",  extract / "extracted_CDS.fa")
    shutil.copy(FIXTURE_DIR / "manifest.tsv",            workdir / "manifest.tsv")
    return SimpleNamespace(
        canonical_gff=extract / "canonical.gff",
        manifest_tsv=workdir / "manifest.tsv",
        extract_dir=extract,
    )


def run_one(paths, model: str, metric_config: dict, output_path: Path, log) -> dict:
    run_core(paths, metric_config=metric_config, output_path=output_path, model=model, log=log)
    with open(output_path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 1:
        raise AssertionError(f"Expected 1 row in {output_path.name}, got {len(rows)}")
    return rows[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--keep-outputs",
        type=Path,
        default=None,
        help="Stage fixtures and outputs in this dir instead of a tempdir (preserved after exit).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Show INFO-level log messages from the plugin.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s %(name)s: %(message)s",
    )
    log = logging.getLogger("nmd_check")

    if args.keep_outputs:
        workdir = args.keep_outputs
        workdir.mkdir(parents=True, exist_ok=True)
        cleanup = None
    else:
        cleanup = tempfile.TemporaryDirectory()
        workdir = Path(cleanup.name)

    try:
        failures: list[str] = []
        n_total = 0
        for strand in ("plus", "minus"):
            print(f"\n--- {strand.upper()} STRAND ---")
            strand_workdir = workdir / strand
            strand_workdir.mkdir(parents=True, exist_ok=True)
            paths = stage_fixtures(strand_workdir, strand)

            for model, cfg, expected in CASES:
                n_total += 1
                label = f"{model}" + (f" [{cfg}]" if cfg else "")
                output_name = model + ("_rule_off" if cfg.get("apply_nmd_rule") is False else "") + ".tsv"
                output_path = strand_workdir / output_name

                row = run_one(paths, model, cfg, output_path, log)
                actual = {k: int(row[k]) for k in expected}

                ok = actual == expected
                marker = "OK  " if ok else "FAIL"
                print(f"  [{marker}] {label:40s}")
                print(f"           expected={expected}")
                print(f"           actual  ={actual}")
                if not ok:
                    diffs = [f"{k}: {expected[k]} vs {actual[k]}" for k in expected if expected[k] != actual[k]]
                    failures.append(f"[{strand}] {label} — " + ", ".join(diffs))

        print()
        if failures:
            print(f"{len(failures)} of {n_total} case(s) FAILED:")
            for f in failures:
                print(f"  - {f}")
            return 1
        print(f"All {n_total} cases passed (across both strands).")
        return 0
    finally:
        if cleanup is not None:
            cleanup.cleanup()


if __name__ == "__main__":
    sys.exit(main())
