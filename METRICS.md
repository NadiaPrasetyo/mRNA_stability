# Metrics

Per-transcript metric tables produced by `bin/01b_metrics.py` and the plugin
modules under `metrics/`. This file is the reference for column names,
units, edge-case handling, and join keys. For how to add a new metric, see
the README.

## Conventions across all metric tables

* **Files**: tab-separated, header on line 1, written to
  `runs/<dataset>/metrics/<plugin>.tsv`.
* **Missing or undefined values**: written as `NA`, consistent with
  `manifest.tsv` and trivially readable by R (`read.delim(..., na.strings="NA")`).
* **Join keys**: every table has `transcript_id` (preserved verbatim from
  `manifest.tsv`, joinable to it directly) and `gene_id` (normalised:
  namespace prefixes like `gene:` / `rna-` stripped, trailing version
  suffix stripped). `gene_id` is suitable for joining across metric tables
  and across datasets if they share an annotation system.
* **Coordinates**: all distances within a transcript are in **spliced
  (mature mRNA) coordinates** unless explicitly labelled as genomic.
  Intron lengths are the one routine exception — necessarily genomic.
* **Direction**: where a metric depends on a 5'→3' reading direction
  (first / last exon, codon-relative distances), it uses transcript
  reading direction, not genomic order. For `-` strand transcripts this
  means the "first exon" is the one with the highest genomic start.
* **U vs T**: T and U are treated equivalently throughout. Output column
  names use the RNA convention (`frac_U` for the combined count).
* **Format**: each plugin chooses long or wide based on its data shape;
  the per-plugin section below documents which.

## Output format choice — long vs wide

Wide is the default (one row per transcript) and is what `junctions` and
`architecture` use. `sequence_basic` is long (one row per transcript per
region) because the same set of composition columns repeats across many
regions; a wide layout would balloon to dozens of region-prefixed columns
most of which are NA for any given analysis.

When joining a long table to a wide one in R:

```
left_join(wide_table, long_table %>% filter(region == "3UTR"),
          by = "transcript_id")
```

or use `pivot_wider(names_from = region)` to flatten if region-prefixed
columns are wanted.

---

## `junctions.tsv`

Exon-junction features per transcript. All distances in spliced
(mature mRNA) coordinates. Wide format.

| column                          | description                                                                                       |
| ------------------------------- | ------------------------------------------------------------------------------------------------- |
| `transcript_id`                 | from `manifest.tsv`, preserved verbatim                                                           |
| `gene_id`                       | normalised                                                                                        |
| `strand`                        | `+` or `-`                                                                                        |
| `n_exons`                       | total exon count                                                                                  |
| `n_5UTR_junctions`              | junctions within 5'UTR features = `len(five_prime_UTR features) - 1`                              |
| `n_CDS_junctions`               | junctions within CDS features = `len(CDS features) - 1`                                           |
| `n_3UTR_junctions`              | junctions within 3'UTR features = `len(three_prime_UTR features) - 1`                             |
| `stop_dist_closest_upstream`    | spliced nt from stop codon to nearest upstream junction; `NA` if none                             |
| `stop_dist_closest_downstream`  | spliced nt from stop codon to nearest downstream junction; `NA` if stop is in last exon           |
| `stop_dist_last_downstream`     | spliced nt to the *last* downstream junction (canonical NMD metric); `NA` if stop is in last exon |
| `start_dist_closest_upstream`   | spliced nt from start codon to nearest upstream junction; `NA` if none                            |
| `start_dist_closest_downstream` | spliced nt from start codon to nearest downstream junction; `NA` if none                          |

**Reference points.** First nt of start codon, last nt of stop codon —
both in transcript reading direction. Falls back to CDS extent if explicit
`start_codon` / `stop_codon` features are absent.

**NMD threshold.** Analyses typically threshold `stop_dist_last_downstream`
at 50 nt. Kept here as a raw distance so threshold sweeps don't need
recomputation.

**Junction sign convention.** A junction located exactly at the reference
position is treated as *upstream* with distance 0. A junction one nt
downstream of the reference is *downstream* with distance 1. This matters
only at the boundary case of a codon coinciding with a splice site.

**Sharp edges:**

* Region junction counts assume the GFF splits UTRs and CDS per-exon (MANE,
  GENCODE, modern Ensembl). Annotations that collapse a UTR into one
  feature spanning introns will report 0 UTR junctions even when the UTR
  spans multiple exons.
* Falling back to CDS extents for the stop codon is off by 3 nt on legacy
  RefSeq annotations that exclude the stop codon from CDS (most modern
  annotations include it). Marginal effect on closest-junction distances.

---

## `architecture.tsv`

Exon and intron length statistics per transcript. Wide format.

| column                 | description                                                                      |
| ---------------------- | -------------------------------------------------------------------------------- |
| `transcript_id`        | from `manifest.tsv`, preserved verbatim                                          |
| `gene_id`              | normalised                                                                       |
| `strand`               | `+` or `-`                                                                       |
| `n_exons`              | total exon count                                                                 |
| `n_internal_exons`     | `max(0, n_exons - 2)`                                                            |
| `first_exon_length`    | nt; transcript reading direction                                                 |
| `last_exon_length`     | nt; transcript reading direction; equals `first_exon_length` when `n_exons == 1` |
| `internal_exon_mean`   | mean length of internal exons (nt); `NA` if `n_internal_exons == 0`              |
| `internal_exon_median` | median length of internal exons (nt); `NA` if `n_internal_exons == 0`            |
| `internal_exon_sd`     | sample SD of internal exon lengths (nt); `NA` if `n_internal_exons < 2`          |
| `intron_mean`          | mean intron length (genomic nt); `NA` if `n_exons < 2`                           |
| `intron_median`        | median intron length (genomic nt); `NA` if `n_exons < 2`                         |
| `intron_sd`            | sample SD of intron lengths (genomic nt); `NA` if `n_exons < 3`                  |

**SD convention.** Sample SD (Bessel-corrected, n-1 denominator) — matches
R's `sd()`. Defined only when the underlying set has at least two values.

**Intron lengths.** Computed as `next_exon.start - this_exon.end - 1` in
genomic order. The resulting set is direction-independent, so set-level
statistics (mean / median / SD) are unaffected by strand. Reported as
genomic, not spliced, nt — they are intron lengths, not transcript distances.

**Single-exon transcripts.** `n_exons = 1`; `first_exon_length` and
`last_exon_length` both equal the exon length; all stat columns are `NA`.

---

## `sequence_basic.tsv`

Length and base composition per transcript per region. Long format — one
row per `(transcript, region)` pair.

| column          | description                                                                                              |
| --------------- | -------------------------------------------------------------------------------------------------------- |
| `transcript_id` | from `manifest.tsv`, preserved verbatim                                                                  |
| `gene_id`       | normalised                                                                                               |
| `region`        | region name (`mRNA`, `CDS`, `5UTR`, `3UTR`, `start_codon_region`, `stop_codon_region`, `tail_region`, …) |
| `length`        | nt                                                                                                       |
| `gc_content`    | `(G+C) / length`; `NA` if `length == 0`                                                                  |
| `frac_A`        | `A / length`; `NA` if `length == 0`                                                                      |
| `frac_C`        | `C / length`; `NA` if `length == 0`                                                                      |
| `frac_G`        | `G / length`; `NA` if `length == 0`                                                                      |
| `frac_U`        | `(T+U) / length`; `NA` if `length == 0`                                                                  |
| `frac_other`    | non-canonical bases (N, R, Y, …) / length; `NA` if `length == 0`                                         |
| `gc_skew`       | `(G-C) / (G+C)`; `NA` if `G+C == 0`                                                                      |
| `at_skew`       | `(A-U) / (A+U)`; `NA` if `A+U == 0`                                                                      |
| `purine_ratio`  | `(A+G) / length`; `NA` if `length == 0`                                                                  |
| `amino_ratio`   | `(A+C) / length`; `NA` if `length == 0`                                                                  |

The five `frac_*` columns sum to 1 by construction (assuming non-zero
length). `gc_content` is computed over total length, not over canonical
bases only — so a heavily-masked sequence gets a lower `gc_content` than
a clean one with the same canonical-base ratio. This makes `frac_other` a
direct read on how much sequence was ambiguous and keeps GC interpretation
transparent.

**Skew conventions.** `gc_skew` and `at_skew` use canonical-only denominators
(standard literature definition; `(G-C)/(G+C)` and `(A-U)/(A+U)`). The two
*ratios* use total length in the denominator (matching the `frac_*` columns).
They are equal to `frac_A + frac_G` and `frac_A + frac_C` respectively,
but exposed as named columns for analysis convenience. Amino bases are
A (NH2 at C6) and C (NH2 at C4); keto bases G and U/T have C=O at the
analogous position — IUPAC code `M` = aMino.

**Region discovery.** The plugin globs `extracted_*.fa` under
`extract_dir`. Region name is the filename's middle (`extracted_<region>.fa`
→ `<region>`). `UTR_pair` is hardcoded as skipped because its 3-line-per-record
ViennaRNA-constraint format is not standard FASTA and would corrupt
composition stats if naively parsed.

**Optional config:**

```
metrics:
  sequence_basic:
    enabled: true
    regions: [mRNA, CDS, 5UTR, 3UTR]    # explicit region list overrides glob
    skip_regions: [tail_region]         # additional skips beyond UTR_pair
```

**Sharp edges:**

* Length-0 records (shouldn't occur in normal extraction output but
  possible if upstream emits them) get `length = 0` and `NA` everywhere
  else.
* T (DNA letter) and U (RNA letter) both count toward `frac_U`, so the
  same plugin runs over DNA-FASTA and RNA-FASTA outputs without change.
* A FASTA header whose first whitespace-token doesn't match any manifest
  transcript_id (after tolerant ID matching) is dropped with a warning.
  Verbatim match should be the norm since extraction writes manifest IDs
  as FASTA headers.

---

## `nmd_fragility_{core,full,window}.tsv`

Density of NMD tripwires in the NMD-competent zones of the CDS:
out-of-frame stop codons (potential ambush stops if a frameshift
occurs upstream) and in-frame codons one single-nucleotide variant
away from a stop codon (potential premature stops if a nonsense
mutation occurs). In-frame fragile codons are split into transition-,
transversion-, and all-SNV-sensitive classes.

Three output files, one per NMD model — identical columns,
identical inputs, different definitions of "competent". Wide format.

| column                               | description                                                                                                                   |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------- |
| `transcript_id`                      | from `manifest.tsv`, preserved verbatim                                                                                       |
| `gene_id`                            | normalised                                                                                                                    |
| `strand`                             | `+` or `-`                                                                                                                    |
| `model`                              | `core`, `full`, or `window`; records which variant produced this row                                                          |
| `cds_length`                         | total extracted CDS length (nt)                                                                                               |
| `zone_length`                        | total length of NMD-competent CDS positions (nt)                                                                              |
| `n_transition_fragile_codons`        | count of in-frame codons one transition away from a stop codon, with the codon's final nucleotide in the NMD-competent zone   |
| `n_transversion_fragile_codons`      | count of in-frame codons one transversion away from a stop codon, with the codon's final nucleotide in the NMD-competent zone |
| `n_snv_fragile_codons`               | count of in-frame codons one transition or transversion away from a stop codon; union of the two fragile classes              |
| `n_alt_stop_codons`                  | count of out-of-frame TAA/TAG/TGA codons whose final nucleotide lies within the NMD-competent zone                            |
| `transition_fragile_codon_density`   | `n_transition_fragile_codons / zone_length`; `NA` if `zone_length == 0`                                                       |
| `transversion_fragile_codon_density` | `n_transversion_fragile_codons / zone_length`; `NA` if `zone_length == 0`                                                     |
| `snv_fragile_codon_density`          | `n_snv_fragile_codons / zone_length`; `NA` if `zone_length == 0`                                                              |
| `alt_stop_codon_density`             | `n_alt_stop_codons / zone_length`; `NA` if `zone_length == 0`                                                                 |
| `transition_fraction_of_snv_fragile` | `n_transition_fragile_codons / n_snv_fragile_codons`; `NA` if `n_snv_fragile_codons == 0`                                     |

**The three models.** The `core` and `full` models threshold the
downstream-junction distance at 50 nt; they differ in *which* junction
matters. The `window` model uses a configurable threshold
(default 50 nt) and restricts scanning to a configurable window.

In simplified terms:

* `core`: scans the core regions of internal exons, excluding the
  final 50 nt before each exon junction.

* `full`: scans the entire contiguous NMD-competent CDS region up to
  the terminal safe zone.

* `window`: scans only a configurable upstream window near the final
  exon-exon junction.

* *Core* (`nmd_fragility_core.tsv`): position `p` is competent iff
  the next downstream junction is more than 50 nt past `p`.
  Equivalent per-exon formulation: within each non-terminal exon,
  positions in `[exon_start, exon_end − 50)` are competent.
  The last 50 nt of each non-terminal exon and the entire terminal
  exon are excluded. Models the assumption that a close downstream
  junction masks the NMD signal from more distal ones.

* *Full* (`nmd_fragility_full.tsv`): position `p` is competent iff
  the last (most-downstream) junction is more than 50 nt past `p`.
  The competent zone is a single contiguous CDS prefix ending at
  `last_junction − 50`. Models the assumption that closer EJCs
  get stripped or occluded by the terminating ribosome anyway, so
  more distal junctions remain and can fire NMD. Matches the
  semantics of `junctions.tsv`'s `stop_dist_last_downstream > 50`.

* *Window* (`nmd_fragility_window.tsv`): evaluates a configurable
  window (default 120 nt) anchored to the last exon-exon junction.
  A configurable threshold (default 50 nt, literature 50-55)
  sets the size of the NMD-immune safe zone, and a boolean toggle
  (`apply_nmd_rule`) determines whether the window respects the
  safe zone or runs right up to the junction. If the available
  upstream CDS is shorter than the window, it scans the available
  sequence and adjusts the denominator accordingly.

**Optional config for the window model:**

```
metrics:
  nmd_fragility_window:
    enabled: true
    window_size: 120          # length of the upstream scanning window (nt)
    nmd_threshold: 50         # NMD-immune safe zone (nt); literature: 50-55
    apply_nmd_rule: true      # set to false to ignore the safe zone
```

The `core` and `full` models diverge on positions in the last 50 nt
of internal exons: `core` excludes them, `full` includes them.
Worked example for a 5-exon CDS with exons of 200 nt each
(junctions at spliced 200, 400, 600, 800; CDS length 1000):

* Core zone:   `[0,150) ∪ [200,350) ∪ [400,550) ∪ [600,750)` = 600 nt
* Full zone:   `[0,750)` = 750 nt
* Window zone (default rules): `[630,750)` = 120 nt

**Fragile codons.** In-frame fragile codons are split by the
substitution class that can convert them to a stop codon.

* Transition-fragile codons: `{CAA, CAG, CGA, TGG}`. These are one
  transition (C↔T or A↔G) away from a stop codon: `CAA→TAA`,
  `CAG→TAG`, `CGA→TGA`, `TGG→TAG`, and `TGG→TGA`.
* Transversion-fragile codons: `{AAA, AAG, AGA, GAA, GAG, GGA, TAC,
  TAT, TCA, TCG, TGC, TGT, TTA, TTG}`. These are one transversion
  away from a stop codon.
* SNV-fragile codons are the union of transition- and
  transversion-fragile codons.

Fragile codons are scanned only in-frame (`i % 3 == 0` relative to the
annotated start codon), where a nonsense mutation would create a
premature termination codon at that codon's position.

**Out-of-frame stops.** Scanned at every CDS position with
`i % 3 != 0`. Positions where `i % 3 == 1` start codons in frame
`+1`, and `i % 3 == 2` start codons in frame `+2`. Models the
possibility that an upstream frameshift would shift the ribosome
onto one of these alternative reading frames, with a downstream
alt-stop then terminating translation.

**Codon-zone convention.** A candidate codon is counted when the final
nucleotide of that codon lies inside the NMD-competent zone. This matches
the `junctions.tsv` reference convention, where stop-codon distances are
measured from the last nucleotide of the stop codon. It also avoids counting
a codon merely because its first base lies inside the zone when the
termination position would fall outside it.

**Density denominators.** All density columns normalise by `zone_length`,
giving burden per NMD-competent nucleotide. `transition_fraction_of_snv_fragile`
is not a density; it is the fraction of all in-frame SNV-fragile codons that
are transition-fragile rather than transversion-only.

**Reference frame.** The annotated start codon defines frame 0.
Its spliced position is taken from `start_codon` features if
present, falling back to the lowest CDS coordinate (`+` strand)
or highest (`-` strand).

**Coverage.** Rows are emitted only for transcripts present in
both `manifest.tsv` and `extracted_CDS.fa`. Non-coding manifest
transcripts do not appear; a `left_join` from `manifest.tsv` in R
preserves them with `NA` rows.

**Sharp edges:**

* Transcripts with no exons, or where the start codon cannot be
  located in spliced coordinates, are skipped (with a counted warning).
* Boundary convention: a junction at exactly position `p` is treated
  as upstream (distance 0), matching `junctions.py`. A junction at
  `p + 1` is downstream with distance 1.
* The threshold is strict (`>`, not `>=`): for the default 50 nt
  threshold, a junction exactly 50 nt downstream of `p` does NOT
  make `p` competent; 51 nt does.

## `uorf.tsv`

Upstream open reading frame (uORF) counts and summary statistics
per transcript. Wide format.

| column                            | description                                                                                   |
| --------------------------------- | --------------------------------------------------------------------------------------------- |
| `transcript_id`                   | from `manifest.tsv`, preserved verbatim                                                       |
| `gene_id`                         | normalised                                                                                    |
| `utr5_length`                     | nt; `NA` if no 5'UTR record                                                                   |
| `n_uorfs`                         | count of classical uORFs (stop within 5'UTR)                                                  |
| `n_overlapping_uorfs`             | count of oORFs (stop at or past main ATG)                                                     |
| `has_uorf`                        | 0/1; 1 iff `n_uorfs + n_overlapping_uorfs > 0`                                                |
| `total_classical_uorf_codons`     | sum of codon counts across classical uORFs                                                    |
| `max_classical_uorf_codons`       | longest classical uORF (codons); `NA` if none                                                 |
| `dist_cap_to_first_uatg`          | nt from 5' end to first uATG (any class); `NA` if none                                        |
| `dist_last_uorf_stop_to_main_atg` | intercistronic nt between last classical uORF's stop and main ATG; `NA` if no classical uORFs |

**Classical vs overlapping.** A uORF is *classical* if its in-frame stop
codon falls entirely within the 5'UTR. It is *overlapping* (oORF) if the
in-frame stop falls at or past the main start codon — typically in the
CDS, terminating at or before the main stop. The two classes have
different functional profiles; oORFs are associated with stronger NMD
effects because their stop codons are deep in the mRNA body.

**Codon-count convention.** Includes the start codon, excludes the stop
(`ATG TAA` = 1 codon; `ATG NNN TAA` = 2 codons).

**Reading frame.** Independent of the main ORF frame. uATGs in any of
the three possible 5'UTR frames are counted, including uATGs that
happen to be in-frame with the main ORF (these almost always produce
oORFs that extend to the main stop codon).

**Start codon.** ATG only by default. The `start_codons` config knob
accepts an explicit list (e.g., `[ATG, CTG, GTG]`) for near-cognate
inclusion.

**Boundary handling.** A uATG must fit entirely within the 5'UTR — i.e.,
all three nt of the start codon are at positions `< utr5_length`. ATGs
straddling the UTR-CDS boundary are not counted (they are essentially
the main ATG with a 1- or 2-nt shift).

**Intercistronic distance.** `dist_last_uorf_stop_to_main_atg` is the
count of nt between the last nt of the classical stop codon and the
first nt of the main ATG — i.e., `utr5_length - (stop_pos + 3)`. A
distance of 0 means the classical uORF stop ends exactly at the main
ATG. Standard "intercistronic distance" definition, predicts
re-initiation efficiency.

**Optional config:**

```
metrics:
  uorf:
    enabled: true
    start_codons: [ATG]      # default; extend with [ATG, CTG, GTG, ACG] for near-cognates
    min_codons: 1            # default; minimum codon count for a uORF to be counted
```

**Sharp edges:**

* Manifest transcript with no 5'UTR FASTA record (typically non-coding,
  or no annotated UTR): NA row, preserving manifest universe. Joins
  cleanly in R.
* 5'UTR present but < 3 nt: zero counts (cannot host a start codon).
* uATG with no in-frame stop within `5'UTR + CDS`: not counted as either
  class. Logged at INFO. Rare in practice; could be extended to scan
  into 3'UTR if it becomes important.
* The scan does not currently weight uORFs by Kozak context strength.
  All uATGs are treated equally regardless of surrounding sequence
  quality. A future extension could add a Kozak-score column.

---

## Joining tables in R

Typical downstream pattern:

```
manifest    <- read.delim("runs/<dataset>/extracted_regions/manifest.tsv", na = "NA")
junctions   <- read.delim("runs/<dataset>/metrics/junctions.tsv", na = "NA")
architect   <- read.delim("runs/<dataset>/metrics/architecture.tsv", na = "NA")
seq_basic   <- read.delim("runs/<dataset>/metrics/sequence_basic.tsv", na = "NA")
nmd_core    <- read.delim("runs/<dataset>/metrics/nmd_fragility_nearest.tsv", na = "NA")
nmd_full     <- read.delim("runs/<dataset>/metrics/nmd_fragility_any.tsv", na = "NA")
nmd_window  <- read.delim("runs/<dataset>/metrics/nmd_fragility_window.tsv", na = "NA")
uorf        <- read.delim("runs/<dataset>/metrics/uorf.tsv", na = "NA")

# Wide table per transcript (long → wide on the basic-composition stats)
library(dplyr); library(tidyr)
seq_wide <- seq_basic %>%
  pivot_wider(id_cols = c(transcript_id, gene_id),
              names_from = region,
              values_from = c(length, gc_content, frac_A, frac_C, frac_G,
                              frac_U, frac_other, gc_skew, at_skew,
                              purine_ratio, amino_ratio))

merged <- manifest %>%
  left_join(junctions,  by = c("transcript_id", "gene_id")) %>%
  left_join(architect,  by = c("transcript_id", "gene_id")) %>%
  left_join(seq_wide,   by = c("transcript_id", "gene_id")) %>%
  left_join(nmd_near,   by = c("transcript_id", "gene_id"),
            suffix = c("", "_core")) %>%
  left_join(nmd_full,    by = c("transcript_id", "gene_id"),
            suffix = c("", "_full")) %>%
  left_join(nmd_window, by = c("transcript_id", "gene_id"),
            suffix = c("", "_window")) %>%
  left_join(uorf,       by = c("transcript_id", "gene_id"))
```

`n_exons` appears in both `junctions.tsv` and `architecture.tsv`. A naive
`left_join` will produce `n_exons.x` / `n_exons.y` columns. Either pick
one source explicitly (`select(-n_exons)` from one table before joining)
or treat them as a consistency check. Similarly, `strand` appears in
`junctions.tsv`, `architecture.tsv`, and all `nmd_fragility_*.tsv`
files — same handling.

**NMD-model handling.** All three `nmd_fragility_*.tsv` files carry an
explicit `model` column. If multiple models are enabled, prefer
loading them as separate tables and joining with suffixes
(e.g., `_core`, `_full`, `_window`) rather than concatenating without
distinguishing, as the rest of the columns share names. Analysts
thresholding `junctions.tsv`'s `stop_dist_last_downstream` at 50 nt
are implicitly applying the full-EJC model to the actual stop codon,
so `nmd_fragility_full` is the natural companion.
