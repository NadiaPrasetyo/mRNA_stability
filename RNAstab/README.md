# RNAstab

A modular R pipeline for analysing mRNA transcript half-life and the sequence,
structural, and translational features that predict it. Designed to run on any
species with an Ensembl-style transcript annotation and the appropriate raw
inputs in place.

Built around the half-life data from Agarwal & Kelley, *The genetic and
biochemical determinants of mRNA degradation rates in mammals*, Genome Biology
23:245 (2022), but agnostic to data provenance — drop in any `halflife.csv`
keyed on `ensembl_gene_id` or `gene_name`.


## Quick start

```r
# 1. Place raw inputs under data/raw/<species>/
# 2. Load the pipeline:
source("R/load_all.R")

# 3. Build (or load from cache) a single species:
human <- build_dataset("human")

# 4. Or build everything you have configured and stack:
all_species <- build_all()      # tibble with a `species` column

# 5. Force a rebuild after changing raw inputs or feature logic:
human <- build_dataset("human", rebuild = TRUE)
```


## Directory layout

```
RNAstab/
├── R/
│   ├── config.R                # paths, species registry, FEATURE_PATTERNS,
│   │                           # CACHE_VERSION
│   ├── load_all.R              # sources every R/ file in dependency order
│   ├── utils/
│   │   ├── normalise.R         # z_score_normalize, min_max_normalize,
│   │   │                       # normalize_numeric
│   │   ├── naming.R            # format_col_name() — canonical → display
│   │   └── feature_groups.R    # fg(), fg_columns(), , select_features(),
|   |                           # resolve_selection()
│   ├── io/
│   │   ├── load_raw.R          # one loader per source file
│   │   └── cache.R             # save_snapshot, load_snapshot, clear_snapshot
│   ├── features/
│   │   ├── mfe_model.R         # calculate_mfe_expected, calculate_mfe_delta
│   │   └── engineer.R          # add_mfe_expected_and_delta, impute_mrna_mfe,
│   │                           # add_junction_density, add_exon_fractions,
│   │                           # engineer_features
│   └── pipeline/
│       ├── assemble.R          # pivot_regional_to_wide, join helpers
│       └── build_dataset.R     # build_dataset(), build_all()
├── analysis/
│   ├── qc/
│   │   ├── mfe_expected_check.R
│   │   └── imputation_check.R
│   ├── correlations/
│   │   ├── scatter_plot.R      # create_scatter_plot()
│   │   └── halflife_correlation.R
│   └── models/
│       └── lasso.R             # run_lasso(), lasso_coefficient_plot()
├── scripts/
│   ├── preprocess_saluki.R     # one-off: HDF5 → .rds
│   ├── build_human.R           # CLI runner for human
│   ├── build_mouse.R           # CLI runner for mouse
│   └── example_analysis.R      # worked end-to-end example
├── data/
│   ├── raw/                    # populated by you
│   │   ├── human/
│   │   ├── mouse/
│   │   └── shared/
│   ├── cache/                  # auto-generated .rds snapshots
│   └── outputs/
│       ├── plots/
│       └── tables/
└── README.md
```


## How it works

### One row per (transcript, species)

Every row of the built dataset is one transcript for one species. A `species`
column makes cross-species analysis trivial:

```r
combined <- bind_rows(build_dataset("human"), build_dataset("mouse"))
combined |> group_by(species) |> summarise(median_hl = median(halflife, na.rm = TRUE))
```

Or just call `build_all()`.

### No species prefix on column names

Columns are named `{metric}_{region}` in lowercase, e.g. `rnafold_zscore_5utr`,
`length_cds`, `gc_3utr`. Half-life is just `halflife`. The species lives in a
column, not in column names. This is what makes the cross-species stack work.

For pretty plot labels, `format_col_name()` (in `R/utils/naming.R`) turns
canonical column names into display strings (`"mfe_zscore_5utr"` →
`"MFE z-score 5' UTR"`). Its sibling `format_group_name()` (in
`R/utils/palettes.R`) does the same for *selection keys* — group, supergroup,
and bundle names (`"rnafold_zscores"` → `"MFE z-score"`, `"nmd_reported"` →
`"NMD (reported)"`) — for plots that label a facet strip or legend by group.
Standalone columns (`cai`, `expression`, …) are labelled as columns by
`format_col_name()`, not by `format_group_name()`.

### Region vocabulary

```
5utr   cds   3utr   mrna   utrpair   last100   start   stop
```

These appear as suffixes on metric columns. `utrpair` is the combined non-coding
regions (5' UTR + 3' UTR)

### Feature groups

Reach for `fg("rnafold_zscores")` instead of typing out region patterns:

```r
df |> select(halflife, fg("rnafold_zscores"), fg("rnal_zscores"))

# What does a group resolve to?
fg_columns(df, "rnafold_zscores")
# [1] "rnafold_zscore_5utr" "rnafold_zscore_cds"  "rnafold_zscore_3utr"
# [4] "rnafold_zscore_mrna" "rnafold_zscore_utrpair"
```

Defined groups (see `R/config.R`):

`lengths`, `gc`, `nmd`, `architecture`, `rnafold_scores`, `rnafold_zscores`,
`mfe_deltas`, `mfe_expected`, `rnal_scores`, `rnal_zscores`, `rnaup`,
`junctions`, `uorfs`, `orfs`, `stopfree`, `skews`, `distances`, `codon_freqs`,
`aa_freqs`, `nuc_ratios`, `probing`.

Add a new one by appending to `FEATURE_PATTERNS` in `R/config.R`.

### Choosing which columns to plot
 
Feature groups answer "what are all the codon columns?" Often you want less
than a whole group — the top few, one named metric, or everything-but-one. That
is *selection intent*, and it lives in a separate layer from the schema so that
narrowing a plot never means editing `FEATURE_PATTERNS`.
 
Three things can name a set of columns:
 
- **A group** — a `FEATURE_PATTERNS` key, e.g. `"codon_freqs"`. The schema.
- **A supergroup** — a coarse family, e.g. `"structure"`, which expands to all
  its member groups. Also schema (see `SUPERGROUPS` in `R/config.R`).
- **A bundle** — a *reusable named selection* you define, e.g. `"nmd_reported"`.
  Intent, not schema (see `GROUP_BUNDLES` in `R/config.R`).
`select_features()` turns any mix of these — plus optional one-off refinements —
into the actual columns present in your dataframe:
 
```r
# A whole supergroup
select_features(df, groups = "structure")
 
# A reusable named subset (defined once in GROUP_BUNDLES)
select_features(df, groups = "nmd_reported")
 
# One-off: keep only two named columns from the nmd group
select_features(df, groups = "nmd",
                pick = list(nmd = c("nmd_snv_fragile_codon_density_mrna",
                                    "nmd_alt_stop_codon_density_mrna")))
 
# One-off: the whole probing group minus one noisy column
select_features(df, groups = "probing",
                drop = list(probing = "gini_nucleoplasm_cds"))
```
 
`pick` is an allow-list (columns added to the group later stay out until you
name them); `drop` removes from the otherwise-whole group (later additions are
included). Use `pick` for a small fixed subset of an open-ended family, `drop`
for "the family minus a couple of members."
 
A **bundle** is just the reusable form of the same idea. Define it once:
 
```r
# in R/config.R
GROUP_BUNDLES <- list(
  nmd_reported = list(
    groups = "nmd",
    pick   = list(nmd = c("nmd_snv_fragile_codon_density_mrna",
                          "nmd_alt_stop_codon_density_mrna"))
  )
)
```
 
…then pass `groups = "nmd_reported"` anywhere a plot accepts `groups`. The
correlation dotplot and the feature/response scatter both understand groups,
supergroups, bundles, and per-call `pick`/`drop`. If you pass both a bundle and
a caller `pick`/`drop` for the same group, the caller wins.
 
No `CACHE_VERSION` bump is ever needed for any of this — it is selection logic,
not feature engineering.

### Caching

After every successful build, the dataset is saved to
`data/cache/{species}_dataset_v{CACHE_VERSION}.rds`. Subsequent calls to
`build_dataset()` return the cache instantly. Three ways to invalidate:

```r
build_dataset("human", rebuild = TRUE)   # one-off rebuild
clear_snapshot("human")                  # delete cache file
# Bump CACHE_VERSION in R/config.R       # invalidates everyone's cache
```

Bump `CACHE_VERSION` whenever feature-engineering logic changes — that's the
mechanism for keeping caches honest after refactors.

### Adding a new species

Three steps:

1. Add an entry to `SPECIES_CONFIG` in `R/config.R` (copy the human or mouse
   block, adjust `dir` and Dani genome patterns).
2. Place raw files under `data/raw/<species_dir>/` matching the filenames the
   loaders expect.
3. `build_dataset("rat")`.

No other code changes needed.

### Adding a new feature

Two patterns depending on what it is.

**A new derived feature** (computed from existing columns): add a function to
`R/features/engineer.R`, then call it in `engineer_features()`. Bump
`CACHE_VERSION`. Done.

**A new raw input source**: add a loader to `R/io/load_raw.R`, then add it to
the appropriate join block in `R/pipeline/build_dataset.R`. If it has a useful
group, append a regex to `FEATURE_PATTERNS`.


## Common analysis recipes

```r
source("R/load_all.R")
df <- build_dataset("human")

# --- Top correlations with half-life ---
source("analysis/correlations/halflife_correlation.R")
p <- halflife_correlation_plot(df, top_n = 20)
print(p)
write.csv(list(plot, table), "data/outputs/tables/halflife_spearman.csv",
          row.names = FALSE)

# --- A custom scatter ---
source("analysis/correlations/scatter_plot.R")
create_scatter_plot(df,
                    x_var = "mfe_zscore_mrna",
                    y_var = "halflife",
                    color_var = "length_mrna",
                    log_color = TRUE,
                    add_density_rings = TRUE)

# --- LASSO model ---
source("analysis/models/lasso.R")
fit <- run_lasso(df)
fit$r_squared
print(lasso_coefficient_plot(fit, top_n = 12))

# --- Custom predictor bundle ---
custom <- list(
  groups  = c("mfe_zscores", "rnal_zscores"),
  columns = c("cai", "translation_efficiency")
)
fit2 <- run_lasso(df, predictor_spec = custom)
```


## Required R packages

Core pipeline: `dplyr`, `tidyr`, `readr`, `purrr`, `stringr`, `tibble`,
`tidyselect`, `rlang`.

Analysis layer: `ggplot2`, `forcats`, `viridis`, `scales`, `glmnet`.

Saluki preprocessing only: `rhdf5` (Bioconductor).


## Behavioural changes from the previous pipeline

1. **mRNA z-score imputation no longer drops transcripts with any missing
   region z-score.** The old `impute_Zscores.R` filtered them out before
   imputing — defeating the purpose. The new `impute_mrna_mfe()` imputes from
   whatever region z-scores are present (length-weighted) and only returns NA
   when no region data exists.

2. **`stopfree_cds = length_cds` placeholder removed.** That assignment was a
   no-op stub; it's gone.

3. **MFE-model parameters are scoped inside `calculate_mfe_expected()`.** Old
   code put `a, b, c, d` as bare globals — collision risk with anything else
   in the session.

4. **Column names are canonical lowercase, no species prefix.** All downstream
   scripts updated to match. `format_col_name()` rewritten for the new schema.


## Known limitations

- Cache invalidation is by integer version, not by hashing the source files.
  If you edit a raw file without bumping `CACHE_VERSION` or passing
  `rebuild = TRUE`, the stale cache wins. (For per-stage caching with automatic
  staleness detection, the {targets} package is the natural next step.)
- The Saluki preprocessing step is one-off and not integrated into the main
  build. If Saluki outputs change, run `scripts/preprocess_saluki.R` again
  before rebuilding.
- Loaders silently skip missing files. Good for incomplete species, bad if you
  expected a file to be picked up. Watch the `skip (missing): …` messages on
  the first build.
