"""Tests for the lifted manifest-indexing helpers in lib.gff.

These tests are committed BEFORE the implementation, to fix the API
shape and resolution-order semantics for the helper lift planned in
`handoff_session2.md`. The module-level skip means this file collects
cleanly under pytest today (alongside the working nmd_fragility tests)
and starts running automatically once `build_manifest_index` and
`resolve_to_manifest_tid` are added to `lib.gff`.

Contract being committed:

    build_manifest_index(manifest_tx_ids: Iterable[str]) -> dict[str, str]
        Index each manifest tid under raw / strip_prefix / normalise_id
        forms; map each form to the original (raw) manifest tid. On
        collisions, first-write-wins (matching the existing
        `_index_keys_for_feature` convention in lib.gff).

    resolve_to_manifest_tid(record_id, index, region=None) -> Optional[str]
        Try, in order:
          1. split_composite_fasta_id(record_id, region) → parsed tid,
             then strip_prefix(tid), then normalise_id(tid)
          2. record_id raw, then strip_prefix(record_id),
             then normalise_id(record_id)
        Return the first hit's mapped value, or None.

Resolution-order rationale: the parsed tid is more specific than the
raw composite header, so it takes precedence. Each form is also
normalised in case the manifest holds a different variant (e.g.
versioned manifest, unversioned FASTA tid, or vice versa).
"""
from __future__ import annotations

import pytest

try:
    from lib.gff import build_manifest_index, resolve_to_manifest_tid
except ImportError:
    pytest.skip(
        "Manifest helpers not yet lifted into lib.gff. "
        "See handoff_session2.md § 'Recommended next move' for the lift plan.",
        allow_module_level=True,
    )


# --- build_manifest_index ---------------------------------------------------

class TestBuildManifestIndex:
    """The index must contain every plausible lookup form for each manifest tid."""

    def test_empty_input(self):
        assert build_manifest_index([]) == {}

    @pytest.mark.parametrize("manifest_id,expected_keys", [
        # Canonical (no prefix, no version) — only one form.
        ("TX001",
         {"TX001"}),
        # Versioned — indexed raw AND version-stripped.
        ("ENST00000259605.11",
         {"ENST00000259605.11", "ENST00000259605"}),
        # Prefixed — indexed raw AND prefix-stripped.
        ("transcript:ENST00000259605",
         {"transcript:ENST00000259605", "ENST00000259605"}),
        # Both prefix and version — all three forms indexed.
        ("transcript:ENST00000259605.11",
         {"transcript:ENST00000259605.11",
          "ENST00000259605.11",
          "ENST00000259605"}),
    ], ids=["canonical", "versioned", "prefixed", "prefixed+versioned"])
    def test_indexes_all_forms(self, manifest_id, expected_keys):
        index = build_manifest_index([manifest_id])
        missing = expected_keys - set(index.keys())
        assert not missing, (
            f"Index missing expected keys for {manifest_id!r}: {missing}. "
            f"Got keys: {sorted(index.keys())}"
        )

    def test_all_keys_map_to_original_tid(self):
        """Each indexed form must map to the *original* manifest string."""
        manifest_id = "transcript:ENST00000259605.11"
        index = build_manifest_index([manifest_id])
        for key in index:
            assert index[key] == manifest_id, (
                f"Key {key!r} maps to {index[key]!r}, expected {manifest_id!r}"
            )

    def test_collision_first_write_wins(self):
        """Two manifest tids that normalise to the same form: first one keeps the slot.

        Matches the existing `_index_keys_for_feature` convention in lib.gff.
        Real annotations don't collide here once the canonical.gff filter has
        trimmed to one transcript per gene, but the behaviour must be defined.
        """
        index = build_manifest_index([
            "ENST00000259605.11",
            "ENST00000259605.10",
        ])
        # The normalised form is shared; first-seen wins.
        assert index["ENST00000259605"] == "ENST00000259605.11"
        # Each raw form still maps to itself.
        assert index["ENST00000259605.11"] == "ENST00000259605.11"
        assert index["ENST00000259605.10"] == "ENST00000259605.10"


# --- resolve_to_manifest_tid ------------------------------------------------

class TestResolveToManifestTid:
    """Resolution must succeed across the (FASTA form × manifest form) cross-product."""

    def test_empty_index_returns_none(self):
        assert resolve_to_manifest_tid("TX001", {}) is None

    def test_unresolvable_returns_none(self):
        index = build_manifest_index(["TX001"])
        assert resolve_to_manifest_tid("TX_NONEXISTENT", index) is None

    def test_canonical_direct_hit(self):
        """FASTA record ID is the bare manifest tid, no composite header."""
        index = build_manifest_index(["TX001"])
        assert resolve_to_manifest_tid("TX001", index) == "TX001"

    def test_composite_with_region_hint(self):
        index = build_manifest_index(["TX001"])
        result = resolve_to_manifest_tid(
            "GENE001_TX001_CDS", index, region="CDS",
        )
        assert result == "TX001"

    def test_composite_multitoken_region(self):
        """Multi-token regions like UTR_pair are only resolvable with the hint."""
        index = build_manifest_index(["TX001"])
        result = resolve_to_manifest_tid(
            "GENE001_TX001_UTR_pair", index, region="UTR_pair",
        )
        assert result == "TX001"

    def test_composite_without_region_hint(self):
        """Without a region hint, the helper falls back to '_'-split."""
        index = build_manifest_index(["TX001"])
        result = resolve_to_manifest_tid("GENE001_TX001_CDS", index)
        assert result == "TX001"

    @pytest.mark.parametrize("manifest_form,fasta_form", [
        # The cross-product verified independently in handoff_session2.md.
        # Every (manifest, FASTA) combination must resolve.
        ("ENST00000259605",                  "ENST00000259605"),
        ("ENST00000259605",                  "ENST00000259605.11"),
        ("ENST00000259605",                  "transcript:ENST00000259605"),
        ("ENST00000259605",                  "transcript:ENST00000259605.11"),
        ("ENST00000259605.11",               "ENST00000259605"),
        ("ENST00000259605.11",               "ENST00000259605.11"),
        ("ENST00000259605.11",               "transcript:ENST00000259605"),
        ("ENST00000259605.11",               "transcript:ENST00000259605.11"),
        ("transcript:ENST00000259605",       "ENST00000259605"),
        ("transcript:ENST00000259605",       "ENST00000259605.11"),
        ("transcript:ENST00000259605",       "transcript:ENST00000259605"),
        ("transcript:ENST00000259605",       "transcript:ENST00000259605.11"),
        ("transcript:ENST00000259605.11",    "ENST00000259605"),
        ("transcript:ENST00000259605.11",    "ENST00000259605.11"),
        ("transcript:ENST00000259605.11",    "transcript:ENST00000259605"),
        ("transcript:ENST00000259605.11",    "transcript:ENST00000259605.11"),
    ])
    def test_cross_product_bare_fasta(self, manifest_form, fasta_form):
        """No composite header — direct cross-product."""
        index = build_manifest_index([manifest_form])
        result = resolve_to_manifest_tid(fasta_form, index)
        assert result == manifest_form, (
            f"Failed: manifest={manifest_form!r}, fasta={fasta_form!r}, "
            f"resolved to {result!r}"
        )

    @pytest.mark.parametrize("manifest_form,fasta_tid", [
        ("ENST00000259605",     "ENST00000259605"),
        ("ENST00000259605",     "ENST00000259605.11"),
        ("ENST00000259605.11",  "ENST00000259605"),
        ("ENST00000259605.11",  "ENST00000259605.11"),
    ])
    def test_cross_product_composite_fasta(self, manifest_form, fasta_tid):
        """Composite FASTA header with region hint — parsed tid takes precedence."""
        record_id = f"GENE001_{fasta_tid}_CDS"
        index = build_manifest_index([manifest_form])
        result = resolve_to_manifest_tid(record_id, index, region="CDS")
        assert result == manifest_form, (
            f"Failed: manifest={manifest_form!r}, record_id={record_id!r}, "
            f"resolved to {result!r}"
        )

    def test_wrong_region_hint_falls_back_gracefully(self):
        """Region hint mismatch should not crash; resolution falls through.

        If a CDS plugin is somehow given a 5UTR-tagged FASTA record, the
        composite parse returns None (region suffix doesn't match), then
        the raw / prefix-stripped / normalised forms of the full composite
        string are tried. None of those match a manifest tid, so the
        result is None — graceful failure, not a crash.
        """
        index = build_manifest_index(["TX001"])
        result = resolve_to_manifest_tid(
            "GENE001_TX001_5UTR", index, region="CDS",
        )
        assert result is None


# --- Integration: build + resolve together ----------------------------------

def test_typical_plugin_usage():
    """End-to-end pattern as it'll appear in plugin code after the lift.

    Gene IDs in composite FASTA headers can't contain underscores
    (split_composite_fasta_id uses partition('_') for the gene/transcript
    split). Transcript IDs and region names may contain underscores.
    """
    manifest_tids = ["ENST00000259605.11", "ENST00000000001.7"]
    index = build_manifest_index(manifest_tids)

    fasta_records = [
        "ENSG00000137075_ENST00000259605_CDS",     # unversioned in FASTA
        "ENSG00000000003_ENST00000000001.7_CDS",   # versioned in FASTA
        "ENSG00000000004_TXNONEXISTENT_CDS",       # not in manifest
    ]
    resolved = [
        resolve_to_manifest_tid(rid, index, region="CDS")
        for rid in fasta_records
    ]
    assert resolved == [
        "ENST00000259605.11",
        "ENST00000000001.7",
        None,
    ]
