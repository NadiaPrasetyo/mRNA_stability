# =============================================================================
# Assembly
# =============================================================================
# Takes the raw-loaded tibbles and produces a single wide-form dataframe:
# one row per transcript, one column per (metric × region) combination,
# plus transcript-level metadata and gene-level attributes.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})


#' Join region-level long dataframes and pivot to wide form
#'
#' Takes a list of long-form tibbles (each with columns transcript_id, region,
#' plus one or more metric columns) and returns a single wide tibble keyed by
#' transcript_id, with columns named `{metric}_{region}`.
#'
#' @param regional_dfs Named list of long-form tibbles (NULLs are dropped).
#' @return A wide-form tibble, or NULL if all inputs were NULL.
pivot_regional_to_wide <- function(regional_dfs) {
  regional_dfs <- Filter(Negate(is.null), regional_dfs)
  if (length(regional_dfs) == 0) return(NULL)

  # Full-join on (transcript_id, region) — any transcript-region pair present
  # in any source is kept, with NAs for missing metrics.
  long <- reduce(regional_dfs, full_join, by = c("transcript_id", "region"))

  # Pivot: one column per metric per region.
  pivot_wider(
    long,
    id_cols    = transcript_id,
    names_from = region,
    values_from = -c(transcript_id, region),
    names_sep  = "_"
  )
}


#' Attach transcript-level features (no region column) to the wide frame
#'
#' Each input tibble is full-joined on `transcript_id`. NULLs are skipped.
join_transcript_level <- function(wide_df, transcript_level_dfs) {
  dfs <- Filter(Negate(is.null), transcript_level_dfs)
  for (df in dfs) {
    if (!"transcript_id" %in% names(df)) next
    wide_df <- left_join(wide_df, df, by = "transcript_id")
  }
  wide_df
}


#' Attach gene-level features (keyed by gene_id)
join_gene_level <- function(wide_df, gene_level_dfs) {
  dfs <- Filter(Negate(is.null), gene_level_dfs)
  for (df in dfs) {
    if (!"gene_id" %in% names(df)) next
    wide_df <- left_join(wide_df, df, by = "gene_id")
  }
  wide_df
}


#' Attach features keyed by gene_name (common for TE, probing data)
join_gene_name_level <- function(wide_df, name_level_dfs) {
  dfs <- Filter(Negate(is.null), name_level_dfs)
  for (df in dfs) {
    if (!"gene_name" %in% names(df)) next
    wide_df <- left_join(wide_df, df, by = "gene_name")
  }
  wide_df
}
