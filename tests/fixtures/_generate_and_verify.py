"""Generate the synthetic CDS FASTA fixture and verify predicted counts.

This script does two things:

1. Builds the 999-nt synthetic CDS sequence used by
   `test_nmd_fragility_synthetic.py` and writes it to
   `synthetic_5x200_CDS.fa`. The sequence is a GCC tandem (no fragile
   codons, no alt-stops in any frame) with targeted inserts at known
   positions.

2. Independently re-implements the competent-zone logic for each of
   the three NMD models and counts the expected number of fragile
   codons and alt-stops. This is a *self-check* on the predictions
   baked into the test file — if this script's predictions ever
   diverge from `test_nmd_fragility_synthetic.py`'s expected values,
   one of them is wrong and the discrepancy must be resolved before
   trusting the test.

Run with `python _generate_and_verify.py`. Exits 0 on success, 1 on
mismatch.
"""
from __future__ import annotations

import sys
from pathlib import Path

# --- Geometry ----------------------------------------------------------------

# Mirrors the 5x200 example documented in METRICS.md.
CDS_LENGTH = 999          # nt; 999 = 333 codons (clean ATG..TAA span)
START_S = 0               # spliced position of start codon
JUNCTIONS = [200, 400, 600, 800]  # spliced positions of exon-exon junctions

# Three NMD models' default parameters.
NMD_THRESHOLD = 50
DISTAL_WINDOW_SIZE = 120

# --- Insert plan -------------------------------------------------------------

# Each entry: (sequence_position, insert_string, comment).
# Positions are 0-indexed spliced coordinates and must be in-frame
# (position % 3 == 0).
#
# CGA insert: contributes 1 fragile codon at i=position (in-frame),
#             0 alt-stops at surrounding out-of-frame positions.
#
# ATAACT insert (6 nt, two in-frame codons "ATA","ACT"): contributes
#             1 alt-stop "TAA" at i=position+1 (out-of-frame),
#             0 fragile codons (ATA and ACT are not in the fragile set).
INSERTS = [
    (99,  "CGA",    "frag@99 — in nearest+any, outside distal_window default"),
    (174, "CGA",    "frag@174 — in any only (close to internal junction)"),
    (660, "CGA",    "frag@660 — in all three default zones"),
    (690, "CGA",    "frag@690 — in all three; survives rule=off shift"),
    (720, "CGA",    "frag@720 — in all three; survives rule=off shift (alt-stop control)"),
    (765, "CGA",    "frag@765 — in none default; only in distal_window rule=off"),
    (297, "ATAACT", "altstop@298 — in nearest+any, outside distal_window default"),
    (675, "ATAACT", "altstop@676 — in all three default zones"),
]

ALT_STOPS = {"TAA", "TAG", "TGA"}
FRAGILE_CODONS = {"CGA", "CAG", "CAA", "TGG"}


# --- Sequence construction ---------------------------------------------------

def build_cds_sequence(length: int = CDS_LENGTH) -> str:
    """Build the synthetic CDS as GCC tandem + ATG/TAA bookends + targeted inserts."""
    assert length % 3 == 0, "CDS length must be a multiple of 3"
    seq = list("GCC" * (length // 3))
    seq[0:3] = list("ATG")          # start codon
    seq[length - 3:length] = list("TAA")  # stop codon

    for pos, insert, _comment in INSERTS:
        assert pos % 3 == 0, f"Insert position {pos} must be in-frame"
        assert pos + len(insert) <= length, f"Insert at {pos} runs past CDS end"
        seq[pos:pos + len(insert)] = list(insert)

    # Verify no two inserts overlap.
    occupied = {}
    for pos, insert, comment in INSERTS:
        for i in range(pos, pos + len(insert)):
            if i in occupied:
                raise AssertionError(
                    f"Insert overlap at position {i}: "
                    f"{occupied[i]!r} vs {comment!r}"
                )
            occupied[i] = comment

    return "".join(seq)


# --- Competent-zone calculators (re-implementations of plugin logic) ---------

def zone_nearest(cds_len: int, start_s: int, junctions: list[int]) -> list[bool]:
    """For each CDS position, competent iff nearest downstream junction > 50 away."""
    is_competent = [False] * cds_len
    ji = 0
    while ji < len(junctions) and junctions[ji] <= start_s:
        ji += 1
    for i in range(cds_len):
        p = start_s + i
        while ji < len(junctions) and junctions[ji] <= p:
            ji += 1
        if ji < len(junctions) and (junctions[ji] - p) > NMD_THRESHOLD:
            is_competent[i] = True
    return is_competent


def zone_any(cds_len: int, start_s: int, junctions: list[int]) -> list[bool]:
    """Competent zone = single prefix ending at last_junction - 50."""
    is_competent = [False] * cds_len
    if not junctions:
        return is_competent
    zone_end = max(0, min(cds_len, junctions[-1] - start_s - NMD_THRESHOLD))
    for i in range(zone_end):
        is_competent[i] = True
    return is_competent


def zone_distal_window(
    cds_len: int,
    start_s: int,
    junctions: list[int],
    window_size: int = DISTAL_WINDOW_SIZE,
    apply_nmd_rule: bool = True,
    nmd_threshold: int = NMD_THRESHOLD,
) -> list[bool]:
    is_competent = [False] * cds_len
    if not junctions:
        return is_competent
    target_j = junctions[-1]
    global_end = target_j - nmd_threshold if apply_nmd_rule else target_j
    global_start = global_end - window_size
    cds_start = max(0, min(cds_len, global_start - start_s))
    cds_end   = max(0, min(cds_len, global_end - start_s))
    for i in range(cds_start, cds_end):
        is_competent[i] = True
    return is_competent


# --- Counter (mirrors plugin's per-transcript counting loop) -----------------

def count_features(cds_seq: str, is_competent: list[bool]) -> tuple[int, int]:
    """Return (n_alt_stops, n_fragile_codons) using the same rules as the plugin."""
    cds_len = len(cds_seq)
    n_alt_stops, n_fragile = 0, 0
    for i in range(cds_len - 2):
        if not is_competent[i]:
            continue
        codon = cds_seq[i:i + 3]
        if i % 3 != 0:
            if codon in ALT_STOPS:
                n_alt_stops += 1
        else:
            if codon in FRAGILE_CODONS:
                n_fragile += 1
    return n_alt_stops, n_fragile


# --- Expected values baked into the test file --------------------------------

# These MUST match the values asserted in test_nmd_fragility_synthetic.py.
EXPECTED = {
    "nearest":                  {"zone_length": 600, "n_fragile": 4, "n_alt_stops": 2},
    "any":                      {"zone_length": 750, "n_fragile": 5, "n_alt_stops": 2},
    "distal_window_default":    {"zone_length": 120, "n_fragile": 3, "n_alt_stops": 1},
    "distal_window_rule_off":   {"zone_length": 120, "n_fragile": 3, "n_alt_stops": 0},
}


# --- Main --------------------------------------------------------------------

def main() -> int:
    here = Path(__file__).parent
    seq = build_cds_sequence()
    assert len(seq) == CDS_LENGTH

    # Write the FASTA.
    fasta_path = here / "synthetic_5x200_CDS.fa"
    header = "GENE001_TX001_CDS"
    # Wrap at 60 columns for readability (standard FASTA convention).
    wrapped = "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))
    fasta_path.write_text(f">{header}\n{wrapped}\n")
    print(f"Wrote {fasta_path} ({len(seq)} nt, header={header!r})")

    # Verify expected counts.
    failures = []

    def check(label, zone, expected):
        zone_length = sum(zone)
        n_alt, n_frag = count_features(seq, zone)
        actual = {"zone_length": zone_length, "n_fragile": n_frag, "n_alt_stops": n_alt}
        ok = actual == expected
        marker = "OK  " if ok else "FAIL"
        print(f"  [{marker}] {label:30s} expected={expected} actual={actual}")
        if not ok:
            failures.append((label, expected, actual))

    print("\nSelf-check of competent-zone calculators + feature counter:")
    check("nearest",
          zone_nearest(CDS_LENGTH, START_S, JUNCTIONS),
          EXPECTED["nearest"])
    check("any",
          zone_any(CDS_LENGTH, START_S, JUNCTIONS),
          EXPECTED["any"])
    check("distal_window (default)",
          zone_distal_window(CDS_LENGTH, START_S, JUNCTIONS),
          EXPECTED["distal_window_default"])
    check("distal_window (rule off)",
          zone_distal_window(CDS_LENGTH, START_S, JUNCTIONS, apply_nmd_rule=False),
          EXPECTED["distal_window_rule_off"])

    # Sanity: sequence-level invariants.
    print("\nSequence invariants:")
    assert seq[0:3] == "ATG", "Start codon must be ATG"
    print(f"  [OK  ] starts with ATG")
    assert seq[CDS_LENGTH - 3:CDS_LENGTH] == "TAA", "Stop codon must be TAA"
    print(f"  [OK  ] ends with TAA")
    assert len(seq) == CDS_LENGTH
    print(f"  [OK  ] length = {CDS_LENGTH}")

    if failures:
        print(f"\n{len(failures)} prediction(s) FAILED. Fixture and test file are out of sync.")
        return 1

    # --- GFF cross-strand check -----------------------------------------
    # Verify that both strand fixtures parse to the same spliced layout.
    # This catches phase / coordinate errors in the hand-written GFFs
    # before they confuse a plugin-level failure.
    print("\nGFF strand-twin check (requires gffutils):")
    try:
        import gffutils
    except ImportError:
        print("  [SKIP] gffutils not installed; can't verify GFFs")
        return 0

    def parse_layout(gff_path: Path):
        db = gffutils.create_db(str(gff_path), ":memory:", force=True,
                                keep_order=True, merge_strategy="merge")
        tx = db["transcript:TX001"]
        exons = list(db.children(tx, featuretype="exon"))
        if tx.strand == "-":
            sorted_e = sorted(exons, key=lambda e: e.start, reverse=True)
        else:
            sorted_e = sorted(exons, key=lambda e: e.start)
        cum = 0
        junctions: list[int] = []
        for e in sorted_e[:-1]:
            cum += (e.end - e.start + 1)
            junctions.append(cum)
        starts = list(db.children(tx, featuretype="start_codon"))
        start_g = (max(s.end for s in starts) if tx.strand == "-"
                   else min(s.start for s in starts))
        # Map start_g to spliced
        spliced_start = 0
        for e in sorted_e:
            if e.start <= start_g <= e.end:
                offset = (e.end - start_g) if tx.strand == "-" else (start_g - e.start)
                start_s = spliced_start + offset
                break
            spliced_start += (e.end - e.start + 1)
        return tx.strand, junctions, start_s

    here = Path(__file__).parent
    plus  = parse_layout(here / "synthetic_5x200.gff3")
    minus = parse_layout(here / "synthetic_5x200_minus.gff3")
    print(f"  plus  strand: junctions={plus[1]}, start_s={plus[2]}")
    print(f"  minus strand: junctions={minus[1]}, start_s={minus[2]}")
    if plus[1:] != minus[1:]:
        print("  [FAIL] Strand twins produce different spliced layouts.")
        return 1
    if plus[1] != [200, 400, 600, 800] or plus[2] != 0:
        print(f"  [FAIL] Plus-strand layout doesn't match documented 5x200 geometry.")
        return 1
    print("  [OK  ] Both strands produce junctions=[200,400,600,800], start_s=0.")

    print("\nAll predictions match. Fixtures are consistent with test expectations.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
