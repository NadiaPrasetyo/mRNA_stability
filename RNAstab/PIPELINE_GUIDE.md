# RNAstab — Pipeline Extension Guide

**Audience:** an AI coding instance (or human contributor) tasked with writing new plotting scripts, refactoring legacy scripts onto this pipeline, or extending the pipeline with new features, loaders, or species.

**Authority:** this document is a rulebook, not a tour. Where it says **MUST** / **MUST NOT** / **DO** / **DO NOT**, treat those as constraints, not suggestions. The README explains *what the pipeline is*. This document explains *how to interface with and extend it without breaking the invariants*.

**Companion reading (in this order):** `README.md` → this document → `R/config.R` → `R/io/load_raw.R` → `R/pipeline/build_dataset.R`. Do not write code before reading all four.

---

## 1. Mental model

### 1.1 Data flow

```
data/raw/<species>/         ──┐
data/raw/shared/            ──┤  load_raw.R (one loader per source file)
                              │   returns NULL on missing files (silent)
                              ▼
                      regional (long)  +  transcript-level (wide)
                              │
                              ▼  pipeline/assemble.R
                       pivot_regional_to_wide()
                              │
                              ▼  join_transcript_level / join_gene_level
                       wide tibble, one row per transcript
                              │
                              ▼  features/engineer.R
                       engineer_features()  (derived columns)
                              │
                              ▼  io/cache.R
                  data/cache/<species>_dataset_v<N>.rds
                              │
                              ▼
                  consumer code: build_dataset(species)
```

### 1.2 Shape invariants

After `build_dataset(species)` returns, the dataframe satisfies **all** of these. Any consumer code MAY rely on them; any extender MUST preserve them:

- **One row per transcript.** `transcript_id` is unique within a species.
- **`species` column is present** and equals the species name on every row.
- **Identifier columns come first:** `species`, `transcript_id`, `gene_id`, `gene_name` (whichever exist).
- **Column names are canonical** (lowercase, snake_case, no species prefix). See §2.
- **Region suffix is the last token** of any regional column: `length_5utr`, `gc_content_cds`, `rnafold_zscore_mrna`.
- **All-NA columns are dropped** by `drop_all_na_columns()` at the end of `engineer_features()`. DO NOT rely on a specific column existing without guarding.
- **`halflife` is the canonical response variable.** No species prefix, no suffix.

### 1.3 Cross-species idiom

```r
combined <- bind_rows(build_dataset("human"), build_dataset("mouse"))
# or
combined <- build_all()
combined |> group_by(species) |> summarise(...)
```

Any plot or analysis that supports multi-species **MUST** facet, group, or filter on the `species` column. DO NOT bake species into column names anywhere.

---

## 2. Canonical schema reference

### 2.1 Region vocabulary

These are the only legal region suffixes. Defined in `REGIONS` in `R/config.R`:

| Suffix     | Meaning                                |
|------------|----------------------------------------|
| `5utr`     | 5' untranslated region                 |
| `cds`      | Coding sequence                        |
| `3utr`     | 3' untranslated region                 |
| `mrna`     | Whole mRNA                             |
| `utrpair`  | 5' UTR + 3' UTR (interaction features) |
| `last100`  | Last 100 nt of CDS                     |
| `start`    | Start codon region                     |
| `stop`     | Stop codon region                      |

DO NOT invent new suffixes. If a feature genuinely needs a new region, add it to `REGIONS` first and document it in this table.

### 2.2 Column-name patterns (canonical → display)

Patterns listed left to right: regex (used by `FEATURE_PATTERNS`) → example column → display string from `format_col_name()`.

| Regex (FEATURE_PATTERNS key)               | Example column                  | Display string                  |
|--------------------------------------------|---------------------------------|---------------------------------|
| `lengths` → `^length_`                     | `length_cds`                    | `Length CDS`                    |
| `gc` → `^gc_`                              | `gc_content_5utr`               | `GC content 5' UTR`             |
| `rnafold_scores`  → `^rnafold_score_`   | `rnafold_score_mrna`            | `MFE mRNA`                      |
| `rnafold_zscores`  → `^rnafold_zscore_` | `rnafold_zscore_5utr`           | `MFE z-score 5' UTR`            |
| `mfe_expected` → `^mfe_expected_`          | `mfe_expected_cds`              | `MFE expected CDS`              |
| `mfe_deltas` → `^mfe_delta_`               | `mfe_delta_3utr`                | `MFE Δ 3' UTR`                  |
| `rnal_scores` → `^rnal_score_`             | `rnal_score_5utr`               | `Local MFE 5' UTR`              |
| `rnal_zscores` → `^rnal_zscore_`           | `rnal_zscore_cds`               | `Local MFE z-score CDS`         |
| `rnaup` → `^rnaup_`                        | `rnaup_score_utrpair`           | `RNAup UTR interactions`        |
| `junctions` → `^junctions_`                | `junctions_count_cds`           | `Junction count CDS`            |
| `distances` → `^eej_dist_`                 | `eej_dist_last_stop`            | `Last EEJ dist. from stop codon`|
| `uorfs` → `^(uorf_\|dist_cap_)`            | `uorf_count_mrna`               | `uORF count mRNA`               |
| `orfs` → `^orf_`                           | `orf_length_5utr`               | `uORF length 5' UTR`            |
| `stopfree` → `^stopfree_`                  | `stopfree_3utr`                 | `Stop-free 3' UTR`              |
| `skews` → `^(gc\|at)_skew_`                | `gc_skew_cds`, `at_skew_cds`    | `GC skew CDS`, `AT skew CDS`    |
| `nmd` → `^nmd_`                            | `nmd_fragile_codon_density_mrna`| `NMD fragile codon density mRNA`|
| `architecture` → `^(intron_\|exon_\|noncoding_)` | `intron_length_mean_mrna` | `Mean intron length mRNA`       |
| `codon_freqs` → `^codon_`                  | `codon_aaa_freq_cds`            | `Codon freq. CDS aaa`           |
| `aa_freqs` → `^aa_`                        | `aa_freq_leu`                   | `AA freq. leu`                  |
| `nuc_ratios` → `^(frac_\|purine_\|amino_)` | `purine_ratio_cds`              | `Purine ratio CDS`              |
| `probing` → `^(shape\|keth)_`              | `shape_score_mrna`              | `icSHAPE score mRNA`            |

Reserved single-token columns: `halflife`, `gene_id`, `gene_name`, `transcript_id`, `species`, `translation_efficiency`, `saluki_prediction`, `prediction_difference`, `cai`, `expression`, `orfexondensity`.

> **v4 note.** Whole-transcript scalar metrics — architecture (`intron_length_mean_mrna`, `exon_length_first_mrna`, …), uORF (`uorf_count_mrna`, …) and NMD fragility (`nmd_fragile_codon_density_mrna`, …) — carry the real `mrna` region suffix. The v3 `transcript` pseudo-region and the `window`/`core`/`full` NMD pseudo-regions have been retired; there is now a single NMD fragility model and "NMD" appears only as a metric-name prefix in the display string.

### 2.3 Display-label round-trip rule

Every column that will appear on a plot axis, legend, or title **MUST** have a `format_col_name()` result that is human-readable (no underscores leaking through, no raw regex tokens, no `NA`).

**Test:** for any new column `x` you add to the pipeline, run `format_col_name("x")` interactively. If the output still contains underscores or looks like a debug string, add a `REPLACEMENTS` entry in `R/utils/naming.R`.

### 2.4 Quick lookup index — every `format_col_name()` rule

Alphabetised reference of every active `REPLACEMENTS` rule in `R/utils/naming.R`. Use this to find the display string for a canonical column without scrolling. **Notation:** `prefix_` = prefix-anchored regex (`^prefix_`); `[exact]` = exact-match regex (`^token$`); `_suffix` = suffix-anchored regex (`_suffix$`). For columns built from a prefix + region (e.g. `length_cds`), the region token is substituted via the rules in §2.1.

#### Prefix and exact-match rules

| Pattern                          | Example column                          | Display result                       |
|----------------------------------|-----------------------------------------|--------------------------------------|
| `aa_freq_`                       | `aa_freq_leu`                           | `AA freq. leu`                       |
| `amino_ratio_`                   | `amino_ratio_cds`                       | `Amino ratio CDS`                    |
| `at_skew_`                       | `at_skew_3utr`                          | `AT skew 3' UTR`                     |
| `cai` [exact]                    | `cai`                                   | `CAI`                                |
| `distance_junction_down_start` [exact] | `distance_junction_down_start`    | `Junction downstream distance (start)` |
| `distance_junction_down_stop` [exact]  | `distance_junction_down_stop`     | `Junction downstream distance (stop)`  |
| `distance_junction_up_start` [exact]   | `distance_junction_up_start`      | `Junction upstream distance (start)`   |
| `distance_junction_up_stop` [exact]    | `distance_junction_up_stop`       | `Junction upstream distance (stop)`    |
| `dist_cap_to_first_uatg_`        | `dist_cap_to_first_uatg_mrna`           | `Dist. cap to first uATG mRNA`       |
| `eej_dist_first_`                | `eej_dist_first_start`                  | `First EEJ dist. from start codon`   |
| `eej_dist_last_`                 | `eej_dist_last_stop`                    | `Last EEJ dist. from stop codon`     |
| `exon_count_internal_`           | `exon_count_internal_mrna`              | `Number of internal exons mRNA`      |
| `exon_length_first_`             | `exon_length_first_mrna`                | `First exon length mRNA`             |
| `exon_length_last_`              | `exon_length_last_mrna`                 | `Last exon length mRNA`              |
| `expression` [exact]             | `expression`                            | `Expression`                         |
| `frac_`                          | `frac_a_cds`                            | `Fraction a CDS`                     |
| `freq_cds_`                      | `freq_cds_aaa`                          | `Codon freq. CDS aaa`                |
| `freq_last10_`                   | `freq_last10_aaa`                       | `Codon freq. last-10 aaa`            |
| `gc_`                            | `gc_cds`                                | `GC CDS`                             |
| `gc_content_`                    | `gc_content_5utr`                       | `GC content 5' UTR`                  |
| `gc_skew_`                       | `gc_skew_3utr`                          | `GC skew 3' UTR`                     |
| `gene_id` [exact]                | `gene_id`                               | `Ensembl gene ID`                    |
| `gene_name` [exact]              | `gene_name`                             | `Gene name`                          |
| `halflife` [exact]               | `halflife`                              | `Half-life`                          |
| `intron_length_mean_`            | `intron_length_mean_mrna`               | `Mean intron length mRNA`            |
| `junctions_count_`               | `junctions_count_cds`                   | `Junction count CDS`                 |
| `junctions_density_`             | `junctions_density_cds`                 | `Junction density CDS`               |
| `junctions_` (fallback)          | `junctions_mrna`                        | `Junctions mRNA`                     |
| `keth_`                          | `keth_score_mrna`                       | `Keth-seq score mRNA`                |
| `length_`                        | `length_cds`                            | `Length CDS`                         |
| `mfe_delta_`                     | `mfe_delta_3utr`                        | `MFE Δ 3' UTR`                       |
| `mfe_expected_`                  | `mfe_expected_cds`                      | `MFE expected CDS`                   |
| `mfe_median_`                    | `mfe_median_5utr`                       | `MFE median 5' UTR`                  |
| `mfe_pval_`                      | `mfe_pval_cds`                          | `MFE p-value CDS`                    |
| `noncoding_length_fraction_`     | `noncoding_length_fraction_mrna`        | `Non-coding length fraction mRNA`    |
| `nmd_alt_stop_count_`            | `nmd_alt_stop_count_mrna`               | `NMD alt-stop count mRNA`            |
| `nmd_alt_stop_density_`          | `nmd_alt_stop_density_mrna`             | `NMD alt-stop density mRNA`          |
| `nmd_fragile_codon_count_`       | `nmd_fragile_codon_count_mrna`          | `NMD fragile codon count mRNA`       |
| `nmd_fragile_codon_density_`     | `nmd_fragile_codon_density_mrna`        | `NMD fragile codon density mRNA`     |
| `nuc_ratio_a`                    | `nuc_ratio_a`                           | `Adenine ratio`                      |
| `nuc_ratio_c`                    | `nuc_ratio_c`                           | `Cytosine ratio`                     |
| `nuc_ratio_g`                    | `nuc_ratio_g`                           | `Guanine ratio`                      |
| `nuc_ratio_u`                    | `nuc_ratio_u`                           | `Uracil ratio`                       |
| `orf_length_`                    | `orf_length_5utr`                       | `uORF length 5' UTR`                 |
| `orf_number_`                    | `orf_number_5utr`                       | `uORF count 5' UTR`                  |
| `orf_percent_length_`            | `orf_percent_length_5utr`               | `uORF % length 5' UTR`               |
| `orfexondensity` [exact]         | `orfexondensity`                        | `ORF-exon density`                   |
| `orfj_density` [exact]           | `orfj_density`                          | `ORF-J density`                      |
| `prediction_difference` [exact]  | `prediction_difference`                 | `Prediction difference`              |
| `purine_ratio_`                  | `purine_ratio_cds`                      | `Purine ratio CDS`                   |
| `rnafold_score_`                 | `rnafold_score_mrna`                    | `MFE mRNA`                           |
| `rnafold_zscore_`                | `rnafold_zscore_5utr`                   | `MFE z-score 5' UTR`                 |
| `rnal_pval_`                     | `rnal_pval_cds`                         | `Local MFE p-value CDS`              |
| `rnal_score_`                    | `rnal_score_5utr`                       | `Local MFE 5' UTR`                   |
| `rnal_zscore_`                   | `rnal_zscore_cds`                       | `Local MFE z-score CDS`              |
| `rnaup_pval_`                    | `rnaup_pval_utrpair`                    | `RNAup p-value UTR interactions`     |
| `rnaup_score_`                   | `rnaup_score_utrpair`                   | `RNAup UTR interactions`             |
| `rnaup_zscore_`                  | `rnaup_zscore_utrpair`                  | `RNAup z-score UTR interactions`     |
| `saluki_prediction` [exact]      | `saluki_prediction`                     | `Saluki prediction`                  |
| `shape_`                         | `shape_score_mrna`                      | `icSHAPE score mRNA`                 |
| `species` [exact]                | `species`                               | `Species`                            |
| `stopfree_`                      | `stopfree_3utr`                         | `Stop-free 3' UTR`                   |
| `transcript_id` [exact]          | `transcript_id`                         | `Transcript ID`                      |
| `translation_efficiency` [exact] | `translation_efficiency`                | `Translation efficiency`             |
| `uorf_count_`                    | `uorf_count_mrna`                       | `uORF count mRNA`                    |
| `uorf_present_`                  | `uorf_present_mrna`                     | `Has uORF mRNA`                      |

#### Region-suffix rules

These fire **last**, after every prefix rule has consumed the leading portion
of the column name. Each matches a single leading separator (space *or*
underscore) plus a region token, end-anchored — `[ _]<region>$`. They combine
with a prefix rule to produce the full display string (e.g. prefix
`nmd_fragile_codon_density_` + region `mrna` → `NMD fragile codon density
mRNA`). The region tokens and their display strings are in the table below.

#### Region tokens

Substituted into prefixed columns by the region-suffix rules above. The eight
tokens are the complete `REGIONS` vocabulary (§2.1).

| Token       | Display string         |
|-------------|------------------------|
| `5utr`      | `5' UTR`               |
| `3utr`      | `3' UTR`               |
| `cds`       | `CDS`                  |
| `mrna`      | `mRNA`                 |
| `utrpair`   | `UTR interactions`     |
| `last100`   | `last 100 nt`          |
| `start`     | `start codon`          |
| `stop`      | `stop codon`           |

**If a column you want to display isn't in this index:** it will fall through to the default behaviour — underscores replaced with spaces, no other transformation. Add a `REPLACEMENTS` entry in `R/utils/naming.R` rather than working around it in the plot.

---

## 3. The public interface — the only things you should call

These are the functions every extender code should rely on. DO NOT reach into internals (`load_*`, `pivot_regional_to_wide`, `engineer_features` directly, etc.) from analysis scripts.

| Function                          | File                                       | Purpose                                          |
|-----------------------------------|--------------------------------------------|--------------------------------------------------|
| `build_dataset(species, rebuild)` | `R/pipeline/build_dataset.R`               | Get the wide dataframe for one species           |
| `build_all(species, rebuild)`     | `R/pipeline/build_dataset.R`               | Stack multiple species, one `species` column     |
| `fg(group)`                       | `R/utils/feature_groups.R`                 | Tidyselect spec for a named group                |
| `fg_columns(df, group)`           | `R/utils/feature_groups.R`                 | Inspect what a group resolves to in this df      |
| `select_features(df, groups, pick, drop)` | `R/utils/feature_groups.R` | Resolve groups/supergroups/bundles + pick/drop → column names |
| `resolve_selection(groups, pick, drop)`   | `R/utils/feature_groups.R` | Normalise a selection → `(groups, pick, drop)` triple         |
| `expand_groups(groups)`                   | `R/utils/feature_groups.R` | Group/supergroup/bundle names → `FEATURE_PATTERNS` keys        |
| `SUPERGROUPS`                             | constant in `R/config.R`   | Coarse group → member-group categorisation                     |
| `GROUP_BUNDLES`                           | constant in `R/config.R`   | Reusable named selections (intent, not schema)                 |
| `format_col_name(x)`              | `R/utils/naming.R`                         | Canonical column name → display string (vectorised) |
| `format_group_name(x, kind)`      | `R/utils/palettes.R`                       | Group / supergroup / bundle key → display string (vectorised) |
| `clear_snapshot(species)`         | `R/io/cache.R`                             | Force next build to rebuild                      |
| `REGIONS`                         | constant in `R/config.R`                   | The legal region suffix vocabulary               |
| `FEATURE_PATTERNS`                | constant in `R/config.R`                   | Group → regex registry                           |
| `SPECIES_CONFIG`                  | constant in `R/config.R`                   | Species registry                                 |
| `OUTPUT_DIR`                      | constant in `R/config.R`                   | `data/outputs` — base path for all outputs       |
| `CACHE_VERSION`                   | constant in `R/config.R`                   | Bump when feature-engineering logic changes      |

Calling `source("R/load_all.R")` once at the top of any script loads all of these.

---

## 4. Hard rules

The following are absolute. Violating them will silently corrupt data, break caching, or make plots unrenderable.

### R1 — Always source `R/load_all.R` first

Every script under `analysis/`, `scripts/`, or anywhere else MUST begin with:

```r
source("R/load_all.R")
```

DO NOT `source()` individual pipeline files. DO NOT redefine pipeline functions locally.

### R2 — Always get data via `build_dataset()`

DO NOT read `.rds` cache files directly. DO NOT call loaders. DO NOT bypass the engineering step. The only legal entry points are `build_dataset(species)` and `build_all()`.

### R3 — Select column groups with `fg()`, never hand-rolled regex

**DO:**
```r
df |> select(halflife, fg("rnafold_zscores"), fg("rnal_zscores"))
```

**DO NOT:**
```r
df |> select(halflife, matches("^rnafold_zscore_"), starts_with("rnal"))
```

If the group you need does not exist in `FEATURE_PATTERNS`, **add it** (see §6.4). DO NOT inline a regex.

### R3a — Select *subsets* through the selection layer, never new schema groups
 
`fg()` (R3) selects a whole group. When you need **less than a whole group** —
the top few members, one or two named columns, or everything-but-one — that is
*selection intent*, and it MUST be expressed through the selection layer, not by
adding a narrower entry to `FEATURE_PATTERNS`.
 
`FEATURE_PATTERNS` and `SUPERGROUPS` are **schema**: one entry per real column
family, each family in exactly one supergroup. DO NOT add subset or alias
entries to them (the deleted `*_some`, `mfe_scores`, `mfe_zscores` keys were
exactly this mistake — a subset masquerading as a family). A subset group
pollutes `expand_groups()`, double-assigns columns in any plot that builds a
column→group map, and forces ad-hoc exclusion lists.
 
Three tools express subset intent:
 
```r
# select_features(): the one entry point for "which columns does this use".
# Accepts groups / supergroups / bundles + per-group pick/drop.
select_features(df, groups = "structure")
select_features(df, groups = "nmd",
                pick = list(nmd = c("nmd_snv_fragile_codon_density_mrna",
                                    "nmd_alt_stop_codon_density_mrna")))
select_features(df, groups = "probing",
                drop = list(probing = "gini_nucleoplasm_cds"))
```
 
- **`pick`** = allow-list. Columns added to the group later stay OUT until
  named. Use for a small fixed subset of an open-ended family.
- **`drop`** = remove from the otherwise-whole group. Later additions are
  INCLUDED. Use for "the family minus a couple of noisy members."
If the subset is **reusable** (you reach for it in more than one place, or it
has a name worth remembering), promote it to a **bundle** (see §6.4a) rather
than retyping `pick`/`drop`. Bundles are intent, defined in `GROUP_BUNDLES`
(`R/config.R`); they are NOT schema and never require a `CACHE_VERSION` bump.
 
**Plot functions that accept `groups` SHOULD also accept `pick` and `drop`** and
resolve all three through `resolve_selection()` + `refine_group_columns()` (see
§6.1, step 3). This keeps selection semantics identical across every plot and is
the only sanctioned way to refine a group inside a plot.

**Labelling a group / supergroup / bundle key.** `format_col_name()` is for
*column* names and produces wrong output on selection keys (`aa_freqs` →
"AA freq. freqs", `gc` → "gc"). When a plot renders a group, supergroup, or
bundle **key** as visible text — a facet strip, a legend, an axis tick keyed by
group — pass it through `format_group_name(key, kind)` instead, where `kind` is
one of `"group"`, `"supergroup"`, `"bundle"`, or `"auto"` (resolve the namespace
with the same supergroup → bundle → group precedence as `resolve_selection()`).
Pass an explicit `kind` when the caller knows it (a plot that always resolves to
group keys passes `kind = "group"`); reserve `"auto"` for mixed-token input.
Display strings live in three maps in `R/utils/palettes.R` —
`FEATURE_GROUP_DISPLAY_NAMES`, `SUPERGROUP_DISPLAY_NAMES`, `BUNDLE_DISPLAY_NAMES`
— seeded where the `tools::toTitleCase` fallback would be wrong. Add a map entry
rather than hardcoding a label in the plot.

**Standalones are columns, not group keys.** `cai`, `translation_efficiency`,
`expression`, and `orfexondensity` are NOT in the group table — they are
labelled by `format_col_name()`. A plot whose group column mixes selection keys
with standalones (e.g. `feature_response_scatter.R`) dispatches per element:
selection keys via `format_group_name()`, everything else via
`format_col_name()`. `group_panel_sweep.R` is the single-group exception and
labels columns only — it does NOT use `format_group_name()`.
 
**Exception — single-group tools.** A tool whose entire premise is "one panel
per schema family" (e.g. `feature_group_panel_sweep()`) operates on raw
`FEATURE_PATTERNS` keys only. It MUST NOT accept supergroups or bundles — those
are multi-group / refined objects that contradict the one-family-per-panel
contract. Such tools may keep a local, documented skip set for high-cardinality
groups (e.g. `DEFAULT_SWEEP_SKIP <- c("codon_freqs", "aa_freqs")`) as
default-view ergonomics; that constant lives in the analysis file, not
`config.R`.

### R4 — Format every plot label through `format_col_name()`

Axis labels, legend titles, plot titles, facet strip labels, summary table column headers, and CSV column labels intended for human eyes **MUST** be passed through `format_col_name()` (or use a scale labeller that calls it: `scale_y_discrete(labels = format_col_name)`).

DO NOT hardcode display strings inside plot functions. If a column needs a better label, edit `REPLACEMENTS` in `R/utils/naming.R`, do not work around it in the plot.

### R5 — Guard every column access

Loaders return `NULL` silently when raw files are missing, and `engineer_features()` drops all-NA columns. A column you expect may not be present. Always guard:

```r
if (!"halflife" %in% names(df)) stop("halflife missing; this species lacks the response variable")
if ("mfe_delta_cds" %in% names(df)) { ... }
```

The single exception is `transcript_id` and `species`, which are pipeline invariants (§1.2).

### R6 — Never bake species into column names

DO NOT create columns like `human_halflife`, `mouse_length_cds`. Species belongs in the `species` column. If you need per-species comparisons, pivot or facet.

### R7 — Bump `CACHE_VERSION` when feature engineering changes

Edit `R/config.R`. Increment `CACHE_VERSION` integer. This invalidates every species' cache on next call to `build_dataset()`.

**When to bump:** any change to `R/features/engineer.R`, any change to a loader's output schema, any change to assembly logic.

**When NOT to bump:** new plotting script, new analysis script, anything under `analysis/`.

For a one-off rebuild without bumping the global version, use `build_dataset(species, rebuild = TRUE)` or `clear_snapshot(species)` first.

### R8 — Outputs go under `OUTPUT_DIR`, never anywhere else

Plots: `data/outputs/plots/<name>.<ext>`. Tables: `data/outputs/tables/<name>.csv`. Use `file.path(OUTPUT_DIR, "plots", ...)`, never a relative path or `getwd()`.

DO NOT create new top-level output directories. DO NOT save to the working directory.

### R9 — Plot functions return `list(plot, table)`

Any function that produces a plot derived from summary statistics (correlations, group means, lasso coefficients, top-N rankings) MUST return:

```r
list(plot = <ggplot>, table = <tibble of the underlying data>)
```

The `table` element is the source of truth that the plot visualises. Consumers can write it to CSV without re-running the computation. Pure scatter plots that visualise raw rows of `df` are an exception — they may return just the ggplot.

### R10 — Exclude derived response variables from response-correlation analyses

When computing correlations against `halflife`, you MUST exclude columns that are derived predictions of half-life — they are circular. The standard exclusion list is:

```r
exclude = c("^saluki_prediction$", "^prediction_difference$")
```

This is the default in `correlate_with_response()`. Any new model-output column added to the pipeline MUST be added to this default exclusion.

### R11 — Long-form region data uses `(transcript_id, region)` keys

Any new long-form regional loader MUST return a tibble with both `transcript_id` and a lowercase `region` column. The `normalise_region()` helper in `R/io/load_raw.R` enforces lowercase. The assembly step pivots on these two keys.

### R12 — Transcript-level wide loaders MUST NOT contain a `gene_id` column

`gene_id` is supplied by `load_transcripts()` only. Other loaders that happen to carry it MUST drop it (`select(-any_of("gene_id"))`) before returning. Failing to do so produces `gene_id.x` / `gene_id.y` suffix collisions on join.

---

## 5. Anti-patterns

Things that look reasonable and will silently break the pipeline or downstream analyses. If you find yourself writing any of these, stop.

| Anti-pattern                                                    | Why it breaks                                       | Correct form                                  |
|-----------------------------------------------------------------|-----------------------------------------------------|-----------------------------------------------|
| `readRDS("data/cache/human_dataset_v2.rds")`                    | Bypasses cache version mgmt, breaks on next bump    | `build_dataset("human")`                      |
| `select(df, matches("^rnafold_zscore_"))`                       | Group definition fragmented across files            | `select(df, fg("rnafold_zscores"))`           |
| `labs(x = "MFE z-score (CDS)")`                                 | Display layer drifts from canonical schema          | `labs(x = format_col_name("rnafold_zscore_cds"))` |
| `df$mfe_delta_cds + df$mfe_delta_3utr`                          | Column may not exist; no guard                      | Guard with `%in% names(df)` or `coalesce`     |
| Saving to `"plots/foo.png"`                                     | Lands in working dir, not `data/outputs/`           | `file.path(OUTPUT_DIR, "plots", "foo.png")`   |
| Hard-coding `c("5utr","cds","3utr")` in a loop                  | Misses `mrna`, `utrpair`, etc.                      | Iterate over `REGIONS`                        |
| Adding `species` filter without `if(species %in% ...)` check    | Silent empty plot when species missing              | `if (!any(df$species == "human")) stop(...)`  |
| Calling `engineer_features()` from an analysis script           | Skips cache, may double-engineer                    | `build_dataset()` does this for you            |
| Mutating column names with `rename_with(toupper)` for display   | Breaks `format_col_name()` round-trip               | Format only at the moment of display          |
| Renaming a column produced by a loader inside `engineer.R`      | Downstream `fg()` patterns break                    | Either rename in the loader, or add new col   |
| Adding `nmd_reported = "^nmd_(snv\|alt)"` to `FEATURE_PATTERNS` | Subset masquerading as a schema family; pollutes `expand_groups`, double-assigns columns | Define it in `GROUP_BUNDLES`, or use per-call `pick`/`drop` |
---

## 6. Extension recipes

Pick the recipe that matches your task. Follow every numbered step.

### 6.1 Adding a new plot or analysis script

**Location:**
- Correlation-style plot → `analysis/correlations/<name>.R`
- Model fit / coefficient plot → `analysis/models/<name>.R`
- QC / diagnostic plot → `analysis/qc/<name>.R`
- Anything else → propose a new subdirectory in your PR description

**Steps:**

1. Create the file. Begin with `source("R/load_all.R")`. Load only the additional packages you need (`ggplot2`, `forcats`, `viridis`, etc.).
2. Define a single primary function. Signature pattern:
   ```r
   <verb>_<noun>_plot <- function(df, ..., formatter = format_col_name) { ... }
   ```
   The first argument MUST be the dataframe; the formatter MUST default to `format_col_name`.
3. Inside the function:
   - Use `fg()` for whole-group selection, or `select_features(df, groups,
     pick, drop)` when the plot exposes group selection to its caller.
   - If the function takes a `groups` argument, also take `pick = list()` and
     `drop = list()`, and resolve all three with `resolve_selection()` +
     `refine_group_columns()` so behaviour matches every other plot (R3a).
   - Filter NA rows on the variables you actually plot.
   - Compute the summary table.
   - Build the ggplot using `format_col_name()` (or `formatter`) for every label.
4. Return `list(plot = p, table = <tibble>)`. (See R9 for the exception.)
5. Append a "top-to-bottom run" block at the bottom of the file, guarded so it only runs when the script is executed directly:
   ```r
   if (sys.nframe() == 0 || identical(environment(), globalenv())) {
     df  <- build_dataset("human")
     out <- <your_function>(df)
     print(out$plot)
     dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
     dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
     ggsave(file.path(OUTPUT_DIR, "plots", "<name>.jpg"),
            plot = out$plot, width = 210, height = 148, units = "mm", dpi = 300)
     write.csv(out$table, file.path(OUTPUT_DIR, "tables", "<name>.csv"), row.names = FALSE)
   }
   ```
6. DO NOT bump `CACHE_VERSION`. Analysis scripts never invalidate the cache.

### 6.2 Refactoring a legacy plotting script

Most legacy scripts will have one or more of these problems. Walk this checklist:

| Symptom                                       | Fix                                                  |
|-----------------------------------------------|------------------------------------------------------|
| Hard-coded column names like `human_halflife` | Strip species prefix → `halflife`. Use `species` column for filtering / faceting. |
| Local `format_col_name_v2()` definition       | Delete it. Use the pipeline's `format_col_name()`. Add missing rules to `REPLACEMENTS`. |
| Reads CSV/RDS from disk directly              | Replace with `build_dataset("<species>")`.           |
| Long hand-rolled `select(matches(...))`       | Replace with `fg()` calls.                           |
| Imputes data inline                           | Move the imputation into `R/features/engineer.R` and bump `CACHE_VERSION`. Plots consume the engineered column. |
| Bare global parameters (`a <- 0.8`, etc.)     | Scope inside the function that uses them. Move thermodynamic constants to `R/features/mfe_model.R`. |
| Writes outputs to working directory           | Use `file.path(OUTPUT_DIR, ...)`.                    |
| Returns just a ggplot from a stats-driven plot| Wrap as `list(plot, table)` per R9.                  |
| Uses `mfe_` legacy alias instead of `rnafold_`| Use the canonical column name in code (`rnafold_score_*`); display will render correctly via `format_col_name()`. |

### 6.3 Adding a new derived feature

A "derived" feature is one computed from columns already in the assembled dataframe.

1. Open `R/features/engineer.R`.
2. Add a function `add_<feature_name>(df) -> df`. The function MUST guard on the presence of every input column it reads (R5). It MUST add columns, not mutate existing ones.
3. Add the call to the pipe inside `engineer_features()`. Order matters: place it after any dependencies and before `drop_all_na_columns()`.
4. If the new column(s) should be selectable as a group, add a regex to `FEATURE_PATTERNS` in `R/config.R`.
5. If the new column(s) have a non-obvious display string, add a rule to `REPLACEMENTS` in `R/utils/naming.R`.
6. **Bump `CACHE_VERSION`** in `R/config.R`.
7. Verify: `build_dataset("human", rebuild = TRUE)` and inspect the new columns. Run `format_col_name("<new_col>")` to confirm the display label.

### 6.4 Adding a new feature group

Simpler. Append to `FEATURE_PATTERNS` in `R/config.R`:

```r
my_new_group = "^my_metric_",
```

No cache bump needed (it's a query helper, not data). Verify with `fg_columns(build_dataset("human"), "my_new_group")`.

### 6.4a Adding a reusable column selection (a bundle)
 
A **bundle** is a named, reusable selection — the right home for what the old
`*_some` groups were trying to be. Use one when a particular subset of a group
(or combination of groups) recurs across plots or modelling.
 
A `GROUP_BUNDLES` entry (in `R/config.R`) is a list with any of:
 
| Field    | Meaning                                                        |
|----------|----------------------------------------------------------------|
| `groups` | character vector of group / supergroup / **other bundle** names |
| `pick`   | named list: group key → columns to KEEP from that group        |
| `drop`   | named list: group key → columns to REMOVE from that group      |
 
A bare character vector is shorthand for `list(groups = <vec>)`.
 
```r
GROUP_BUNDLES <- list(
  # "the two reported NMD columns" — a fixed allow-list against an open group
  nmd_reported = list(
    groups = "nmd",
    pick   = list(nmd = c("nmd_snv_fragile_codon_density_mrna",
                          "nmd_alt_stop_codon_density_mrna"))
  ),
  # "the four canonical length columns"
  lengths_core = list(
    groups = "lengths",
    pick   = list(lengths = c("length_5utr", "length_cds",
                              "length_3utr", "length_mrna"))
  )
)
```
 
**Steps:**
 
1. Add the entry to `GROUP_BUNDLES` in `R/config.R`.
2. Use it anywhere `groups` is accepted: `select_features(df, "nmd_reported")`,
   or a plot call `groups = c("structure", "nmd_reported")`.
3. **DO NOT bump `CACHE_VERSION`.** A bundle is a query helper, not data.
4. Verify: `select_features(build_dataset("human"), "nmd_reported")` returns the
   expected columns.
**Merge rule (know this).** If a caller passes `pick`/`drop` AND names a bundle
carrying `pick`/`drop` for the *same group key*, the **caller wins** for that
key; the bundle's entry applies only where the caller is silent. Resolution is
always pick-then-drop, so a caller `drop` trims whatever a bundle `pick`
produced.
 
**`pick` vs `drop`, restated for bundles:** prefer `pick` when the keep-set is
small and fixed and you do NOT want future group additions to appear (the NMD
case). Prefer `drop` when you want the whole family minus a few members,
including anything added later.
 
**Resolution precedence** inside `expand_groups()` / `resolve_selection()` is
**supergroup → bundle → plain group**, first match wins. DO NOT give a bundle
the same name as an existing group or supergroup.

### 6.5 Adding a new raw input source

1. Add a loader to `R/io/load_raw.R`. Required behaviours:
   - Takes a single `species` argument.
   - Calls `read_if_exists(species_path(species, "<filename>"))`.
   - Returns `NULL` if missing (the helper does this).
   - Returns canonical lowercase snake_case column names.
   - If regional/long-form: includes `transcript_id` and lowercase `region` (use `normalise_region()`).
   - If wide-form: includes `transcript_id` and drops `gene_id` (R12).
2. Wire it into `R/pipeline/build_dataset.R`:
   - Long-form → add to the `regional` list.
   - Transcript-level wide-form → add to the `transcript_level` list.
   - Gene-level → add a `load_*` call and an explicit `left_join` block with the right key.
3. If columns form a meaningful group, add a regex to `FEATURE_PATTERNS` (§6.4).
4. Add display rules to `REPLACEMENTS` for any non-obvious column names.
5. Bump `CACHE_VERSION`.

### 6.6 Adding a new species

1. Append a block to `SPECIES_CONFIG` in `R/config.R`. Copy the `human` block, change `dir`, adjust `dani.shape_pattern` / `dani.keth_pattern` / `dani.extra_cols` to match the genome suffixes in the shared probing file.
2. Place raw files under `data/raw/<dir>/` matching the names the loaders expect (`grep -h "species_path" R/io/load_raw.R` lists every expected filename).
3. Add a runner: `scripts/build_<species>.R`. Copy `scripts/build_human.R` verbatim, change one string.
4. `Rscript scripts/build_<species>.R`. Watch the `skip (missing): ...` messages — they tell you which files are absent.

No code changes to any other file should be needed. If they are, you've found a leak — fix the leak rather than patching around it.

---

## 7. Worked examples

Two worked examples covering different plot shapes. The first is correlation-shaped (every column in a group vs the response); the second is distribution-shaped (response distribution stratified by bins of a predictor). Together they exercise most of the public interface and rules.

### 7.1 Correlation panel — every column in a feature group vs `halflife`

**Task:** for a chosen feature group, produce a faceted scatter panel showing every column in the group plotted against `halflife`, with Spearman correlation annotated on each panel.

This example exercises: `build_dataset()`, `fg_columns()`, vectorised `format_col_name()`, the `list(plot, table)` return contract, the runner-block pattern, and `OUTPUT_DIR`.

File: `analysis/correlations/feature_group_panel.R`

```r
# =============================================================================
# Feature-group correlation panel
# =============================================================================
# For a chosen feature group, plot every member column against halflife as a
# small-multiples scatter panel, annotated with Spearman correlation.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_group_panel.R")
#   df  <- build_dataset("human")
#   out <- feature_group_panel(df, group = "rnafold_zscores")
#   print(out$plot)
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})


#' Faceted scatter panel of every column in a feature group vs halflife
#'
#' @param df      A dataframe from build_dataset().
#' @param group   Character. A key of FEATURE_PATTERNS.
#' @param response Character. Response column (default "halflife").
#' @param formatter Function. Display formatter (default format_col_name).
#' @return list(plot, table). `table` has columns variable, n, spearman, p_value.
#' @export
feature_group_panel <- function(df,
                                group,
                                response  = "halflife",
                                formatter = format_col_name) {

  # --- Rule R5: guard every column access -----------------------------------
  if (!response %in% names(df)) {
    stop("response '", response, "' not in df — wrong species or missing data")
  }

  # --- Rule R3: use fg_columns to enumerate ---------------------------------
  cols <- fg_columns(df, group)
  if (length(cols) == 0) {
    stop("Feature group '", group, "' resolved to zero columns in this df")
  }

  # --- Long-format for faceting ---------------------------------------------
  long <- df |>
    dplyr::select(dplyr::all_of(c(response, cols))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(cols),
                        names_to = "variable",
                        values_to = "value") |>
    dplyr::filter(!is.na(value), !is.na(.data[[response]]))

  # --- Summary table (Rule R9: this is what the plot visualises) ------------
  summary_tbl <- long |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      n        = dplyr::n(),
      spearman = suppressWarnings(
        cor(value, .data[[response]], method = "spearman",
            use = "pairwise.complete.obs")),
      p_value  = suppressWarnings(
        cor.test(value, .data[[response]], method = "spearman",
                 exact = FALSE)$p.value),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      p_adj_bh   = p.adjust(p_value, method = "BH"),
      annotation = sprintf("ρ = %.2f\nq = %.2g", spearman, p_adj_bh)
    )

  # --- Rule R4: every label goes through the formatter ----------------------
  facet_labels <- setNames(formatter(summary_tbl$variable),
                           summary_tbl$variable)

  p <- ggplot(long, aes(x = value, y = .data[[response]])) +
    geom_point(alpha = 0.3, size = 0.6, shape = 16) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "#4B0082", linewidth = 0.6) +
    geom_text(
      data = summary_tbl,
      aes(x = -Inf, y = Inf, label = annotation),
      hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE
    ) +
    facet_wrap(~ variable, scales = "free_x",
               labeller = labeller(variable = facet_labels)) +
    labs(
      title    = paste0("Feature group: ", group, " vs ", formatter(response)),
      subtitle = sprintf("%d columns, Spearman with BH-adjusted q",
                         nrow(summary_tbl)),
      x        = NULL,
      y        = formatter(response)
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text    = element_text(size = 9),
      panel.grid.minor = element_blank()
    )

  # --- Rule R9: return both halves ------------------------------------------
  list(plot = p, table = summary_tbl)
}


# --- Top-to-bottom run (Rule §6.1 step 5) ------------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  df  <- build_dataset("human")
  out <- feature_group_panel(df, group = "rnafold_zscores")
  print(out$plot)

  # Rule R8: outputs under OUTPUT_DIR
  dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

  ggsave(file.path(OUTPUT_DIR, "plots", "feature_group_panel_rnafold_zscores.jpg"),
         plot = out$plot, width = 297, height = 210, units = "mm", dpi = 300)
  write.csv(out$table,
            file.path(OUTPUT_DIR, "tables", "feature_group_panel_rnafold_zscores.csv"),
            row.names = FALSE)
}
```

Things to notice in the example:

- **Every rule reference is inline as a comment** so the reader can trace which rule motivates which line.
- **`fg_columns()` is used to enumerate**, not `fg()` — because we need the names as strings for the `summary_tbl$variable` column. `fg()` is for tidyselect, `fg_columns()` is for string vectors.
- **The summary table is what gets written to CSV** — that's the data source of truth, not the plot.
- **`labeller = labeller(variable = facet_labels)`** is the idiomatic way to apply `format_col_name()` to facet strips.

### 7.2 Distribution by predictor bins — `halflife` stratified by a chosen column

**Task:** plot the distribution of `halflife` across bins of a chosen predictor. Auto-detect whether the predictor is numeric (→ quantile bins) or categorical (→ native levels). Run a Kruskal-Wallis test for differences across bins. Facet by `species` automatically if the input contains more than one.

This example demonstrates a different shape from §7.1: a single predictor (not a group), distributions rather than scatter, categorical-vs-numeric input handling, and multi-species facetting via the `species` column (R6).

File: `analysis/correlations/halflife_by_predictor.R`

```r
# =============================================================================
# Halflife distribution stratified by a predictor
# =============================================================================
# Plot the distribution of halflife across bins of a chosen predictor:
#   - numeric predictor → quantile bins (default 4 bins = quartiles)
#   - categorical / logical / factor predictor → its native levels
# Tests whether halflife differs across bins via Kruskal-Wallis.
# Multi-species inputs are faceted by species (Rule R6).
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/halflife_by_predictor.R")
#
#   # Categorical:
#   df  <- build_dataset("human")
#   out <- halflife_by_predictor(df, predictor = "uorf_present_mrna")
#
#   # Numeric → quartiles:
#   out <- halflife_by_predictor(df, predictor = "mfe_delta_cds", n_bins = 4)
#
#   # Multi-species facet:
#   all <- build_all()
#   out <- halflife_by_predictor(all, predictor = "uorf_present_mrna")
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})


#' Halflife distribution stratified by a predictor (numeric → quantile bins,
#' categorical → native levels). Multi-species inputs are faceted by species.
#'
#' @param df         A dataframe from build_dataset() or build_all().
#' @param predictor  Character. Column name of the predictor.
#' @param n_bins     Integer. Number of quantile bins for numeric predictors.
#'                   Ignored for categorical predictors.
#' @param response   Character. Response column (default "halflife").
#' @param formatter  Function. Display formatter (default format_col_name).
#' @return list(plot, table). table columns: species (if present), bin, n,
#'   median, q25, q75, kw_chi2, kw_p, kw_p_adj_bh.
#' @export
halflife_by_predictor <- function(df,
                                  predictor,
                                  n_bins    = 4,
                                  response  = "halflife",
                                  formatter = format_col_name) {

  # --- Rule R5: guard every column access ----------------------------------
  if (!response  %in% names(df)) stop("response '",  response,  "' not in df")
  if (!predictor %in% names(df)) stop("predictor '", predictor, "' not in df")

  has_species <- "species" %in% names(df) &&
                 length(unique(df$species)) > 1

  # --- Auto-detect predictor type & bin if numeric -------------------------
  is_numeric_pred <- is.numeric(df[[predictor]]) &&
                     length(unique(stats::na.omit(df[[predictor]]))) >= n_bins

  if (is_numeric_pred) {
    qs <- stats::quantile(df[[predictor]],
                          probs = seq(0, 1, length.out = n_bins + 1),
                          na.rm = TRUE)
    if (length(unique(qs)) < n_bins + 1) {
      message("Quantile breaks not unique for '", predictor,
              "' — collapsing to ", length(unique(qs)) - 1, " bins")
      qs <- unique(qs)
    }
    df$.bin <- cut(df[[predictor]], breaks = qs, include.lowest = TRUE,
                   labels = paste0("Q", seq_len(length(qs) - 1)))
  } else {
    df$.bin <- as.factor(df[[predictor]])
  }

  # --- Drop missing --------------------------------------------------------
  df <- df |> dplyr::filter(!is.na(.bin), !is.na(.data[[response]]))
  if (nrow(df) == 0) stop("No non-NA rows for predictor / response combination")

  # --- Per-bin summary (Rule R9: this is the data the plot visualises) -----
  group_vars <- if (has_species) c("species", ".bin") else ".bin"
  summary_tbl <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n      = dplyr::n(),
      median = stats::median(.data[[response]], na.rm = TRUE),
      q25    = stats::quantile(.data[[response]], 0.25, na.rm = TRUE),
      q75    = stats::quantile(.data[[response]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(bin = .bin)

  # --- Kruskal-Wallis: per species if faceting, else overall ---------------
  kw_one <- function(sub) {
    res <- suppressWarnings(stats::kruskal.test(sub[[response]] ~ sub$.bin))
    tibble::tibble(kw_chi2 = unname(res$statistic), kw_p = res$p.value)
  }
  if (has_species) {
    kw_tbl <- df |>
      dplyr::group_by(species) |>
      dplyr::group_modify(~ kw_one(.x)) |>
      dplyr::ungroup() |>
      dplyr::mutate(kw_p_adj_bh = stats::p.adjust(kw_p, method = "BH"))
    summary_tbl <- dplyr::left_join(summary_tbl, kw_tbl, by = "species")
  } else {
    kw <- kw_one(df)
    summary_tbl$kw_chi2     <- kw$kw_chi2
    summary_tbl$kw_p        <- kw$kw_p
    summary_tbl$kw_p_adj_bh <- kw$kw_p  # only one test, no correction
  }

  # --- Plot (Rule R4: every label through formatter) -----------------------
  subtitle <- if (has_species) {
    "Kruskal-Wallis per species (BH-adjusted q in facet annotation)"
  } else {
    sprintf("Kruskal-Wallis: chi^2 = %.2f, p = %s",
            summary_tbl$kw_chi2[1],
            format.pval(summary_tbl$kw_p[1], digits = 2, eps = 0.001))
  }

  title_suffix <- if (is_numeric_pred) " (quantile bins)" else ""

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .bin, y = .data[[response]])) +
    ggplot2::geom_violin(fill = "#c8d8e4", alpha = 0.7, scale = "width") +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA,
                          fill = "white", alpha = 0.9) +
    ggplot2::labs(
      title    = paste0(formatter(response), " by ",
                        formatter(predictor), title_suffix),
      subtitle = subtitle,
      x        = formatter(predictor),
      y        = formatter(response)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title    = ggplot2::element_text(size = 13),
      axis.text     = ggplot2::element_text(size = 11),
      panel.grid.minor = ggplot2::element_blank()
    )

  # --- Multi-species facet (Rule R6) ---------------------------------------
  if (has_species) {
    p <- p + ggplot2::facet_wrap(~ species)
    kw_labels <- summary_tbl |>
      dplyr::distinct(species, kw_chi2, kw_p_adj_bh) |>
      dplyr::mutate(label = sprintf("chi^2 = %.2f\nq = %.2g",
                                    kw_chi2, kw_p_adj_bh))
    p <- p + ggplot2::geom_text(
      data = kw_labels,
      ggplot2::aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE
    )
  }

  list(plot = p, table = summary_tbl)
}


# --- Top-to-bottom run (Rule §6.1 step 5) ------------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  df <- build_dataset("human")

  # Example 1: categorical predictor
  out_cat <- halflife_by_predictor(df, predictor = "uorf_present_mrna")
  print(out_cat$plot)

  # Example 2: numeric predictor → quartiles
  out_num <- halflife_by_predictor(df, predictor = "mfe_delta_cds", n_bins = 4)
  print(out_num$plot)

  dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

  ggsave(file.path(OUTPUT_DIR, "plots", "halflife_by_uorf_present_mrna.jpg"),
         plot = out_cat$plot, width = 210, height = 148, units = "mm", dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "plots", "halflife_by_mfe_delta_cds_quartiles.jpg"),
         plot = out_num$plot, width = 210, height = 148, units = "mm", dpi = 300)

  write.csv(out_cat$table,
            file.path(OUTPUT_DIR, "tables", "halflife_by_uorf_present_mrna.csv"),
            row.names = FALSE)
  write.csv(out_num$table,
            file.path(OUTPUT_DIR, "tables", "halflife_by_mfe_delta_cds_quartiles.csv"),
            row.names = FALSE)
}
```

Things to notice that differ from §7.1:

- **Predictor type auto-detection** — `is.numeric()` plus a uniqueness check decides whether to bin or factorise. Plot functions that accept "any column" should always branch on the column's actual type, not assume.
- **Quantile binning with collapse fallback** — when the predictor has heavy ties (e.g. lots of zeros), quantile breaks can be non-unique. The example detects this and reduces the bin count rather than failing.
- **Hypothesis test stored in the table, not just the subtitle.** Per R9, the table must contain the data the plot visualises. The Kruskal-Wallis result is per-facet, so it's joined onto the summary table as columns that repeat across bins within a facet. A downstream consumer reading the CSV can recover the test without re-running it.
- **`species` faceting is conditional** — only applied when the dataframe actually contains multiple species. Single-species inputs get a flat plot. This is how R6 plays out in practice: code branches on the *content* of the `species` column, never on column names.
- **Two `ggsave()` calls in the runner block.** A script may produce multiple outputs; just give each a distinctive filename under `OUTPUT_DIR/plots/`.
- **`stats::` and `ggplot2::` prefixes** on namespace-ambiguous functions (`quantile`, `median`, `aes`) — safe practice in any script that may run in environments where many packages are attached.

---

## 8. Pre-submit checklist

Before considering any extension complete, run through this list. Each item maps to a numbered rule.

- [ ] Script starts with `source("R/load_all.R")` (R1)
- [ ] Data acquired via `build_dataset()` only (R2)
- [ ] All column-group selection uses `fg()` or `fg_columns()` (R3)
- [ ] All plot labels go through `format_col_name()` (R4)
- [ ] Every column access is guarded (R5)
- [ ] No species names in column names (R6)
- [ ] `CACHE_VERSION` bumped iff pipeline logic changed (R7)
- [ ] All outputs under `OUTPUT_DIR` (R8)
- [ ] Stats-driven plot returns `list(plot, table)` (R9)
- [ ] Response-correlation analyses exclude derived predictions (R10)
- [ ] Long-form loaders use `(transcript_id, region)` (R11)
- [ ] Wide-form loaders drop `gene_id` (R12)
- [ ] Any new column has a working `format_col_name()` result (no leftover underscores)
- [ ] Any new feature group has a `FEATURE_PATTERNS` entry
- [ ] Column *subsets* use `pick`/`drop` or a `GROUP_BUNDLES` bundle — never a new `FEATURE_PATTERNS` subset entry (R3a)
- [ ] Plot functions accepting `groups` also accept `pick`/`drop` and resolve via `resolve_selection()` (R3a)
- [ ] Any new bundle name does not collide with a group or supergroup name (§6.4a)
- [ ] Script runs cleanly from a fresh R session via `Rscript <file>`
- [ ] Outputs are written to the correct subdirectory under `data/outputs/`
- [ ] Display-name spot check: `format_col_name(c(<your new cols>))` returns clean strings

---

## 9. Reference: where things live

```
R/                                  pipeline core (DO NOT scatter analysis here)
├── config.R                        constants, paths, FEATURE_PATTERNS, CACHE_VERSION
├── load_all.R                      sources everything in dependency order
├── utils/                          pure helpers, no pipeline state
│   ├── normalise.R                 z_score_normalize, min_max_normalize
│   ├── naming.R                    format_col_name, REPLACEMENTS
│   └── feature_groups.R            fg, fg_columns
├── io/
│   ├── load_raw.R                  one function per raw source file
│   └── cache.R                     save/load/clear snapshot
├── features/
│   ├── mfe_model.R                 thermodynamic constants + math
│   └── engineer.R                  derived features, engineer_features()
└── pipeline/
    ├── assemble.R                  long→wide pivot, join helpers
    └── build_dataset.R             build_dataset, build_all

analysis/                           consumer code (this is where new plots go)
├── qc/                             diagnostics, sanity checks
├── correlations/                   anything correlation- or scatter-shaped
└── models/                         fitted models, coefficient plots

scripts/                            CLI runners (one per build, one-off jobs)
data/
├── raw/                            inputs you place
├── cache/                          .rds snapshots, auto-generated
└── outputs/                        all generated artifacts
    ├── plots/
    └── tables/
```

---

## 10. Known gotchas (read once, remember forever)

- **`build_dataset()` is cached by integer version, not by file hash.** Edit a raw input file without bumping `CACHE_VERSION` or passing `rebuild = TRUE`, and the stale cache silently wins.
- **Loaders skip missing files silently.** Watch the `skip (missing): ...` messages on first build — they're the only signal that an expected input wasn't found.
- **Numeric colour variables with fewer than 10 unique values are auto-discretised** in `create_scatter_plot()`. Don't be surprised when an integer `species_index`-style column shows up as a factor.
- **`add_mfe_expected_and_delta()` keys on `gc_content_{region}` and `length_{region}`.** If your loader produces `gc_{region}` instead, no expected/delta columns are computed. (Current state: `sequence_basic.tsv` produces `gc_content_*`.)
- **`engineer_features()` drops all-NA columns at the end.** A column whose loader produced only `NA` for this species will simply not appear in the output. This is by design but easy to forget when debugging "where did my column go".
- **The `species` column is added before `engineer_features()` is called.** Any engineering step that operates row-wise has access to `species` if it needs species-specific behaviour. Use this sparingly — most logic should be species-agnostic.
- **v4 retired the pseudo-region tokens.** Whole-transcript scalars (architecture, uORF, NMD) now end in the real `mrna` suffix; the v3 `transcript` token and the `window`/`core`/`full` NMD tokens no longer exist. There is one NMD fragility model — `load_nmd_fragility(species)` takes no `type` argument and reads `nmd_fragility.tsv`. Any analysis script written against the v3 names (`*_transcript`, `nmd_*_window`, …) will silently select nothing — `fg()` returns an empty set rather than erroring.
- **Subset feature groups do not belong in `FEATURE_PATTERNS`.** The deleted `*_some` / `mfe_scores` / `mfe_zscores` keys were subsets and aliases living in the schema layer; they double-assigned columns and forced ad-hoc exclusion lists. Reusable subsets are `GROUP_BUNDLES` bundles; one-off subsets are per-call `pick`/`drop`. Neither needs a `CACHE_VERSION` bump.
     
---

## 11. When in doubt

Read the source file first (it's ~50–200 lines per module, all heavily commented). Run `fg_columns(df, "<group>")` and `format_col_name("<column>")` to inspect what the schema actually produces in the current dataframe. The pipeline is small enough to read end-to-end in under an hour — do that before writing anything that fights against it.
