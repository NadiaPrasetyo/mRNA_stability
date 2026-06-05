# =============================================================================
# Pipeline configuration
# =============================================================================
# Central configuration for the RNA half-life analysis pipeline.
# Edit values here to add species, change paths, or register new feature groups.
# =============================================================================

# --- Paths -------------------------------------------------------------------

DATA_ROOT  <- "data"
RAW_DIR    <- file.path(DATA_ROOT, "raw")
SHARED_DIR <- file.path(RAW_DIR, "shared")
CACHE_DIR  <- file.path(DATA_ROOT, "cache")
OUTPUT_DIR <- file.path(DATA_ROOT, "outputs")

# Bump this integer when feature-engineering logic changes so stale caches
# are regenerated instead of silently reused.
CACHE_VERSION <- 7L


# --- Region vocabulary -------------------------------------------------------
# Canonical internal region names are lowercase. Display names for plots are
# handled by `format_col_name()` in R/utils/naming.R.
#
# Regions are an ordered 5' -> 3' traversal of the transcript. `mrna` is the
# whole-mRNA member of regional families (length_mrna = length_5utr+cds+3utr)
# and, as of v4, is also the suffix used for genuinely whole-transcript
# scalar metrics (architecture, uORF, probing, NMD) — the v3 `transcript`
# pseudo-region and the `window`/`core`/`full` NMD pseudo-regions have been
# retired. There are no pseudo-regions: every token here is a real region.

REGIONS <- c("5utr", "cds", "3utr", "mrna", "utrpair",
             "last100", "start", "stop")

#' Raw → canonical region-token aliases.
#'
#' Some upstream pipelines emit verbose / suffixed region names. This is the
#' single mapping `normalise_region()` consults to bring them in line with
#' REGIONS. Add an entry when a new variant appears upstream; do NOT scatter
#' equivalent renames through individual loaders.
REGION_ALIASES <- c(
  tail_region        = "last100",
  start_codon_region = "start",
  stop_codon_region  = "stop"
)


# --- Species registry --------------------------------------------------------
# Add a new species by appending an entry. Each entry defines where to find
# its raw data and which columns to pull from shared files.

SPECIES_CONFIG <- list(
  
  human = list(
    dir        = "human",
    saluki_rds = "saluki_predictions.rds"    # relative to species dir; NULL if none
  ),
  
  mouse = list(
    dir        = "mouse",
    saluki_rds = "saluki_predictions.rds"
  )
)


# --- Feature groups ----------------------------------------------------------
# Regex patterns defining named groups of columns. Use downstream with
# `fg()` / `fg_columns()`. Add a new group by appending a named entry.
#

FEATURE_PATTERNS <- list(
  lengths        = "^length_",
  gc             = "^gc_",
  nmd            = "^nmd_",
  architecture   = "^(intron_|exon_|noncoding_)", 
  mfe_scores     = "^rnafold_score_",
  mfe_zscores    = "^rnafold_zscore_",
  rnafold_scores = "^rnafold_score_",
  rnafold_zscores = "^rnafold_zscore_",
  rnafold_per_nt = "^rnafold_per_nt_",
  mfe_deltas     = "^mfe_delta_",
  mfe_expected   = "^mfe_expected_",
  rnalfold_scores    = "^rnalfold_score_", 
  rnalfold_zscores   = "^rnalfold_zscore_",
  rnaup          = "^rnaup_",
  junctions      = "^(junctions_|eej_dist_)", 
  uorfs          = "^(uorf_|dist_cap_)", 
  orfs           = "^orf_",
  stopfree       = "^stopfree_",
  skews          = "^(gc|at)_skew_",
  distances      = "^eej_dist_",
  codon_freqs    = "^codon_",
  aa_freqs       = "^aa_",
  nuc_ratios     = "^(frac_|purine_|amino_)",
  probing        = "^gini_"
)


# --- Feature supergroups -----------------------------------------------------
# Coarse-grained categorisation of FEATURE_PATTERNS keys for plotting, palette
# assignment, and any analysis that wants to colour or facet by category
# family. Each FEATURE_PATTERNS key SHOULD belong to exactly one supergroup.
# Standalone reserved columns (cai, translation_efficiency, expression,
# orfexondensity) are NOT listed here — they are handled per-plot.
#
# Update this when adding new groups to FEATURE_PATTERNS.

SUPERGROUPS <- list(
  structure  = c("rnafold_scores", "rnafold_zscores", "rnafold_per_nt",
                 "mfe_deltas", "mfe_expected",
                 "rnalfold_scores", "rnalfold_zscores",
                 "rnaup", "probing"),

  intrinsic  = c("lengths", "gc", "stopfree", "skews", "codon_freqs", "aa_freqs", "nuc_ratios"),

  splicing   = c("junctions", "distances", "architecture", "nmd"),

  regulatory = c("uorfs", "orfs")
)

GROUP_BUNDLES <- list(
  nmd_reported = list(
    groups = "nmd",
    pick = list(nmd = c("nmd_snv_fragile_codon_density_mrna",
                        "nmd_alt_stop_codon_density_mrna"))
  ),
  lengths_core = list(
    groups = "lengths",
    pick = list(lengths = c("length_5utr", "length_cds",
                            "length_3utr", "length_mrna"))
  )
)

# --- Helpers -----------------------------------------------------------------

#' Return the absolute path to a raw file for a given species.
#' @param species Character, one of names(SPECIES_CONFIG).
#' @param filename Character, filename within the species folder.
species_path <- function(species, filename) {
  stopifnot(species %in% names(SPECIES_CONFIG))
  file.path(RAW_DIR, SPECIES_CONFIG[[species]]$dir, filename)
}

#' Return the absolute path to a shared raw file.
shared_path <- function(filename) file.path(SHARED_DIR, filename)

#' Return the cache file path for a species.
cache_path <- function(species) {
  file.path(CACHE_DIR, sprintf("%s_dataset_v%d.rds", species, CACHE_VERSION))
}

#' Prefix payload columns with a tool or other prefix, leaving keys untouched.
#' Retained for backward compatibility; new code should prefer affix_payload().
prefix_payload <- function(df, prefix, keys = c("transcript_id", "region")) {
  dplyr::rename_with(df, ~ paste0(prefix, "_", .x), -dplyr::all_of(keys))
}

#' Affix payload columns with a prefix and/or suffix, leaving key columns
#' untouched. Generalises prefix_payload — used where a loader needs to push
#' a token to the END of the column name (e.g. a trailing region suffix)
#' rather than the front.
#'
#' @param df     A dataframe / tibble.
#' @param prefix Character prepended to every non-key column (default "").
#' @param suffix Character appended to every non-key column (default "").
#' @param keys   Character vector of key columns to leave untouched. Keys not
#'   present in `df` are ignored (unlike prefix_payload, which errors).
#' @return df with non-key columns renamed.
affix_payload <- function(df, prefix = "", suffix = "",
                          keys = c("transcript_id", "region")) {
  keep <- intersect(keys, names(df))
  dplyr::rename_with(df, ~ paste0(prefix, .x, suffix), -dplyr::all_of(keep))
}


# --- Supergroup helpers ------------------------------------------------------

#' Reverse lookup: which supergroup does a FEATURE_PATTERNS key belong to?
#' @param group Character vector of FEATURE_PATTERNS keys.
#' @return Character vector of supergroup names; NA_character_ for unknown keys.
#' @examples
#' supergroup_of("rnafold_zscores")   # "structure"
#' supergroup_of(c("gc", "nmd"))      # c("intrinsic", "splicing")
supergroup_of <- function(group) {
  vapply(group, function(g) {
    hits <- names(SUPERGROUPS)[vapply(
      SUPERGROUPS, function(members) g %in% members, logical(1)
    )]
    if (length(hits) == 0) NA_character_ else hits[1]
  }, character(1), USE.NAMES = FALSE)
}


