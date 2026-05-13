"""Sanity check for the three NMD fragility models on a hand-built fixture.

Three assertion layers, each tied to a specific class of potential bug:

* **Layer 1 — zone-length geometry.** For the 5-exon × 200-nt example
  documented in METRICS.md, the three models produce different
  `nmd_zone_length` values (600 / 750 / 120). These are pure functions
  of the spliced-coordinate machinery and the model's competent-zone
  definition; they are independent of sequence content. A bug in
  threshold inequality, junction iteration, or window arithmetic
  shows up here first.

* **Layer 2 — feature counting.** Given a known CDS sequence built
  with targeted CGA and ATAACT inserts at known positions, each model
  produces predictable counts of `n_fragile_codons` and `n_alt_stops`.
  A bug in frame handling, iteration boundary (`range(cds_len - 2)`),
  or codon-set membership shows up here. The three models disagree on
  several inserts (notably i=174, which lands in `any` only), so a
  divergence between expected and actual immediately localises the
  fault.

* **Layer 3 — configurable toggle.** Re-running the distal_window
  model with `apply_nmd_rule: False` shifts the competent zone from
  `[630, 750)` to `[680, 800)`, changing which inserts are inside
  it. This verifies that `metric_config` is being threaded all the
  way to `_competent_zone_distal_window`.

The fixture and expected values are kept in sync by
`fixtures/_generate_and_verify.py`, which re-implements the
competent-zone logic independently and recomputes the predictions.
Run it after editing the fixture to catch drift before pytest does.
"""
from __future__ import annotations

import csv
import logging
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

# Resolved via conftest.py's sys.path injection.
from metrics._nmd_fragility_core import run as run_core  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures"


# --- Fixtures ---------------------------------------------------------------

@pytest.fixture
def synthetic_paths(tmp_path):
    """Copy synthetic fixtures into tmp_path and return a paths-like namespace.

    Mirrors the directory layout the real pipeline gives the plugin:
    `paths.canonical_gff`, `paths.manifest_tsv`, `paths.extract_dir`
    (containing `extracted_CDS.fa` and a gffutils DB cache).
    """
    extract_dir = tmp_path / "extract"
    extract_dir.mkdir()

    shutil.copy(FIXTURE_DIR / "synthetic_5x200.gff3",     extract_dir / "canonical.gff")
    shutil.copy(FIXTURE_DIR / "synthetic_5x200_CDS.fa",   extract_dir / "extracted_CDS.fa")
    shutil.copy(FIXTURE_DIR / "manifest.tsv",             tmp_path / "manifest.tsv")

    return SimpleNamespace(
        canonical_gff=extract_dir / "canonical.gff",
        manifest_tsv=tmp_path / "manifest.tsv",
        extract_dir=extract_dir,
    )


@pytest.fixture
def log():
    logger = logging.getLogger("test_nmd_fragility_synthetic")
    logger.setLevel(logging.INFO)
    return logger


def _read_single_row(tsv_path: Path) -> dict:
    """Read a one-transcript TSV and return its row as a dict.

    Asserts that there is exactly one data row — useful for catching
    plugin bugs that silently emit zero or multiple rows for a
    single-transcript fixture.
    """
    with open(tsv_path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert len(rows) == 1, f"Expected 1 data row in {tsv_path.name}, got {len(rows)}"
    return rows[0]


# --- Layer 1: zone-length geometry ------------------------------------------

@pytest.mark.parametrize(
    "model,expected_zone_length",
    [
        ("nearest",       600),
        ("any",           750),
        ("distal_window", 120),
    ],
)
def test_zone_length(synthetic_paths, tmp_path, log, model, expected_zone_length):
    """Layer 1: nmd_zone_length matches the geometry documented in METRICS.md."""
    output = tmp_path / f"{model}.tsv"
    run_core(synthetic_paths, metric_config={}, output_path=output, model=model, log=log)
    row = _read_single_row(output)
    actual = int(row["nmd_zone_length"])
    assert actual == expected_zone_length, (
        f"{model}: nmd_zone_length mismatch. "
        f"Expected {expected_zone_length} (per METRICS.md's 5x200 example), got {actual}."
    )


# --- Layer 2: feature counting within the competent zone --------------------

@pytest.mark.parametrize(
    "model,expected_fragile,expected_alt_stops",
    [
        # See fixtures/_generate_and_verify.py for the position-by-position
        # derivation of these counts.
        ("nearest",       4, 2),
        ("any",           5, 2),
        ("distal_window", 3, 1),
    ],
)
def test_feature_counts(
    synthetic_paths, tmp_path, log,
    model, expected_fragile, expected_alt_stops,
):
    """Layer 2: feature counts match positions placed in/out of each model's zone."""
    output = tmp_path / f"{model}.tsv"
    run_core(synthetic_paths, metric_config={}, output_path=output, model=model, log=log)
    row = _read_single_row(output)

    actual_fragile = int(row["n_fragile_codons"])
    actual_alt_stops = int(row["n_alt_stops"])

    assert actual_fragile == expected_fragile, (
        f"{model}: n_fragile_codons mismatch. Expected {expected_fragile}, got {actual_fragile}. "
        f"This suggests a bug in in-frame codon iteration or the competent-zone boundary "
        f"for this model."
    )
    assert actual_alt_stops == expected_alt_stops, (
        f"{model}: n_alt_stops mismatch. Expected {expected_alt_stops}, got {actual_alt_stops}. "
        f"This suggests a bug in out-of-frame scanning or codon-set membership."
    )


# --- Layer 3: apply_nmd_rule toggle -----------------------------------------

def test_distal_window_apply_nmd_rule_false(synthetic_paths, tmp_path, log):
    """Layer 3: apply_nmd_rule=False shifts the window and re-scopes counts.

    Default window (rule on):  [630, 750), contains frag@660, frag@690,
        frag@720 and altstop@676 → n_fragile=3, n_alt_stops=1.
    Rule-off window:           [680, 800), contains frag@690, frag@720,
        frag@765 and no alt-stops → n_fragile=3, n_alt_stops=0.

    Zone length is invariant under the toggle (window size unchanged,
    only its position shifts). The discriminating signal is the
    alt-stop count dropping 1→0 and the *identity* of the fragile
    codons shifting (660 out, 765 in) — a bug that fails to act on
    the config flag would leave both at the default values.
    """
    output = tmp_path / "distal_window_rule_off.tsv"
    run_core(
        synthetic_paths,
        metric_config={"apply_nmd_rule": False},
        output_path=output,
        model="distal_window",
        log=log,
    )
    row = _read_single_row(output)
    assert int(row["nmd_zone_length"])  == 120
    assert int(row["n_fragile_codons"]) == 3
    assert int(row["n_alt_stops"])      == 0


# --- Smoke test: fixture is well-formed -------------------------------------

def test_fixture_well_formed(synthetic_paths):
    """Cheap sanity check that the fixture files exist and look right.

    Catches accidental deletion / corruption of fixture files before
    the test reports a confusing plugin-side failure.
    """
    # FASTA: composite header, 999 nt of sequence, starts ATG, ends TAA.
    fa = (synthetic_paths.extract_dir / "extracted_CDS.fa").read_text().splitlines()
    assert fa[0] == ">GENE001_TX001_CDS", f"Unexpected FASTA header: {fa[0]!r}"
    seq = "".join(line.strip() for line in fa[1:])
    assert len(seq) == 999, f"Expected 999 nt CDS, got {len(seq)}"
    assert seq[:3] == "ATG", f"CDS should start ATG, starts with {seq[:3]!r}"
    assert seq[-3:] == "TAA", f"CDS should end TAA, ends with {seq[-3:]!r}"

    # GFF: parseable, has all the expected feature types.
    import gffutils
    db = gffutils.create_db(
        str(synthetic_paths.canonical_gff), ":memory:",
        force=True, keep_order=True, merge_strategy="merge",
    )
    assert len(list(db.features_of_type("gene")))        == 1
    assert len(list(db.features_of_type("mRNA")))        == 1
    assert len(list(db.features_of_type("exon")))        == 5
    assert len(list(db.features_of_type("CDS")))         == 5
    assert len(list(db.features_of_type("start_codon"))) == 1
