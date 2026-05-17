# METRICS.md extract (just the relevent parts for nmd_fragility)

# Metrics

Per-transcript metric tables produced by `bin/01b_metrics.py` and the plugin
modules under `metrics/`. This file is the reference for column names,
units, edge-case handling, and join keys. For how to add a new metric, see
the README.

## Conventions across all metric tables

- **Files**: tab-separated, header on line 1, written to
  `runs/<dataset>/metrics/<plugin>.tsv`.
- **Missing or undefined values**: written as `NA`, consistent with
  `manifest.tsv` and trivially readable by R (`read.delim(..., na.strings="NA")`).
- **Join keys**: every table has `transcript_id` (preserved verbatim from
  `manifest.tsv`, joinable to it directly) and `gene_id` (normalised:
  namespace prefixes like `gene:` / `rna-` stripped, trailing version
  suffix stripped). `gene_id` is suitable for joining across metric tables
  and across datasets if they share an annotation system.
- **Coordinates**: all distances within a transcript are in **spliced
  (mature mRNA) coordinates** unless explicitly labelled as genomic.
  Intron lengths are the one routine exception — necessarily genomic.
- **Direction**: where a metric depends on a 5'→3' reading direction
  (first / last exon, codon-relative distances), it uses transcript
  reading direction, not genomic order. For `-` strand transcripts this
  means the "first exon" is the one with the highest genomic start.
- **U vs T**: T and U are treated equivalently throughout. Output column
  names use the RNA convention (`frac_U` for the combined count).
- **Format**: each plugin chooses long or wide based on its data shape;
  the per-plugin section below documents which.

## Output format choice — long vs wide

Wide is the default (one row per transcript) and is what `junctions` and
`architecture` use. `sequence_basic` is long (one row per transcript per
region) because the same set of composition columns repeats across many
regions; a wide layout would balloon to dozens of region-prefixed columns
most of which are NA for any given analysis.

When joining a long table to a wide one in R:

    left_join(wide_table, long_table %>% filter(region == "3UTR"),
              by = "transcript_id")

or use `pivot_wider(names_from = region)` to flatten if region-prefixed
columns are wanted.

---


## `nmd_fragility_{nearest,any,distal_window}.tsv`

Density of NMD tripwires in the NMD-competent zones of the CDS:
out-of-frame stop codons (potential ambush stops if a frameshift
occurs upstream) and in-frame codons one transition away from a
stop codon (potential premature stops if a nonsense mutation
occurs).

Three output files, one per NMD model — identical columns,
identical inputs, different definitions of "competent". Wide format.

| column | description |
|---|---|
| `transcript_id` | from `manifest.tsv`, preserved verbatim |
| `gene_id` | normalised |
| `strand` | `+` or `-` |
| `model` | `nearest`, `any`, or `distal_window`; records which variant produced this row |
| `cds_length` | total extracted CDS length (nt) |
| `zone_length` | total length of NMD-competent CDS positions (nt) |
| `n_alt_stops` | count of out-of-frame TAA/TAG/TGA codons starting within the NMD-competent zone |
| `n_fragile_codons` | count of in-frame CGA/CAG/CAA/TGG codons starting within the NMD-competent zone |
| `alt_stop_density` | `n_alt_stops / zone_length`; `NA` if `zone_length == 0` |
| `fragile_codon_density` | `n_fragile_codons / zone_length`; `NA` if `zone_length == 0` |

**The three models.** The `nearest` and `any` models threshold the
downstream-junction distance at 50 nt; they differ in *which* junction
matters. The `distal_window` model uses a configurable threshold
(default 50 nt) and restricts scanning to a configurable window.

- *Nearest-EJC* (`nmd_fragility_nearest.tsv`): position `p` is
  competent iff the next downstream junction is more than 50 nt
  past `p`. Equivalent per-exon formulation: within each non-terminal
  exon, positions in `[exon_start, exon_end − 50)` are competent.
  The last 50 nt of each non-terminal exon and the entire terminal
  exon are excluded. Models the assumption that a close downstream
  junction masks the NMD signal from more distal ones.

- *Any-EJC* (`nmd_fragility_any.tsv`): position `p` is competent
  iff the last (most-downstream) junction is more than 50 nt past
  `p`. The competent zone is a single contiguous CDS prefix ending
  at `last_junction − 50`. Models the assumption that closer EJCs
  get stripped or occluded by the terminating ribosome anyway, so
  more distal junctions remain and can fire NMD. Matches the
  semantics of `junctions.tsv`'s `stop_dist_last_downstream > 50`.

- *Distal-Window* (`nmd_fragility_distal_window.tsv`): evaluates a
  configurable window (default 120 nt) anchored to the last exon-exon
  junction. A configurable threshold (default 50 nt, literature 50-55)
  sets the size of the NMD-immune safe zone, and a boolean toggle
  (`apply_nmd_rule`) determines whether the window respects the
  safe zone or runs right up to the junction. If the available
  upstream CDS is shorter than the window, it scans the available
  sequence and adjusts the denominator accordingly.

**Optional config for the distal window:**

    metrics:
      nmd_fragility_distal_window:
        enabled: true
        window_size: 120          # length of the proximal scanning window (nt)
        nmd_threshold: 50         # NMD-immune safe zone (nt); literature: 50-55
        apply_nmd_rule: true      # set to false to ignore the safe zone

The nearest and any models diverge on positions in the last 50 nt of
internal exons: nearest excludes them, any includes them. Worked example
for a 5-exon CDS with exons of 200 nt each (junctions at spliced 200,
400, 600, 800; CDS length 1000):

- Nearest zone: `[0,150) ∪ [200,350) ∪ [400,550) ∪ [600,750)` = 600 nt
- Any zone:     `[0, 750)` = 750 nt
- Distal window (default rules): `[630, 750)` = 120 nt

**Fragile codons.** The set `{CGA, CAG, CAA, TGG}` comprises every
codon that is exactly one transition (C↔T or A↔G) away from a stop
codon: `CGA→TGA`, `CAG→TAG`, `CAA→TAA`, `TGG→TAG`, `TGG→TGA`.
Scanned only in-frame (`i % 3 == 0` relative to the annotated start
codon), where a nonsense mutation would create a premature
termination codon at that codon's position.

**Out-of-frame stops.** Scanned at every CDS position with
`i % 3 != 0`. Positions where `i % 3 == 1` start codons in frame
`+1`, and `i % 3 == 2` start codons in frame `+2`. Models the
possibility that an upstream frameshift would shift the ribosome
onto one of these alternative reading frames, with a downstream
alt-stop then terminating translation.

**Reference frame.** The annotated start codon defines frame 0.
Its spliced position is taken from `start_codon` features if
present, falling back to the lowest CDS coordinate (`+` strand)
or highest (`-` strand).

**Coverage.** Rows are emitted only for transcripts present in
both `manifest.tsv` and `extracted_CDS.fa`. Non-coding manifest
transcripts do not appear; a `left_join` from `manifest.tsv` in R
preserves them with `NA` rows.

**Sharp edges:**
- Transcripts with no exons, or where the start codon cannot be
  located in spliced coordinates, are skipped (with a counted warning).
- Boundary convention: a junction at exactly position `p` is treated
  as upstream (distance 0), matching `junctions.py`. A junction at
  `p + 1` is downstream with distance 1.
- The threshold is strict (`>`, not `>=`): for the default 50 nt
  threshold, a junction exactly 50 nt downstream of `p` does NOT
  make `p` competent; 51 nt does.

---
