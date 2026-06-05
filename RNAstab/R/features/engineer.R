# =============================================================================
# Feature engineering
# =============================================================================
# All derived features live here. Each function takes a wide-form dataframe
# (one row per transcript) and returns the same shape with new columns added.
#
# Operations are composed by `engineer_features()` at the bottom, which is
# the public entry point called from the pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# --- MFE expected / delta for every region ----------------------------------

#' Add `mfe_expected_{region}` and `mfe_delta_{region}` columns
#'
#' Only applied to regions that have all three of: `rnafold_score_{region}`,
#' `gc_content_{region}`, `length_{region}`.
add_mfe_expected_and_delta <- function(df, regions = REGIONS) {
  for (r in regions) {
    score_col <- paste0("rnafold_score_", r)
    gc_col    <- paste0("gc_content_", r)  # Updated naming
    len_col   <- paste0("length_", r)
    exp_col   <- paste0("mfe_expected_", r)
    del_col   <- paste0("mfe_delta_", r)

    if (all(c(score_col, gc_col, len_col) %in% names(df))) {
      df[[exp_col]] <- calculate_mfe_expected(df[[gc_col]], df[[len_col]])
      df[[del_col]] <- calculate_mfe_delta(df[[score_col]], df[[exp_col]])
    }
  }
  df
}


# --- mRNA-level imputation from region scores -------------------------------

#' Impute mRNA-level MFE z-score and score from region-level data
#'
#' When `rnafold_zscore_mrna` is missing but region-level z-scores are available,
#' compute a length-weighted average over the available regions.
#'
#' Fix vs. old pipeline: the old code filtered out any transcript missing a
#' region z-score before imputing — which defeats the point of imputation.
#' This version only requires ≥1 region to be present. Transcripts with no
#' region z-scores at all get NA, as before.
#'
#' Adds columns:
#'   - rnafold_zscore_mrna_imputed : length-weighted avg of region z-scores
#'   - rnafold_score_mrna_imputed  : sum of region scores (pure additive)
#' Then coalesces the real mRNA columns with the imputed ones:
#'   - rnafold_score_mrna, rnafold_zscore_mrna now filled for transcripts that had
#'     region-level but not mRNA-level data.
#' Finally recomputes mfe_expected_mrna / mfe_delta_mrna on the coalesced
#' score, in case any rows were newly filled.
impute_mrna_mfe <- function(df) {
  required_z   <- paste0("rnafold_zscore_", c("5utr", "cds", "3utr"))
  required_sc  <- paste0("rnafold_score_",  c("5utr", "cds", "3utr"))
  required_len <- paste0("length_",         c("5utr", "cds", "3utr"))
  
  have_z   <- all(required_z   %in% names(df))
  have_sc  <- all(required_sc  %in% names(df))
  have_len <- all(required_len %in% names(df))
  
  if (!(have_z && have_len)) {
    message("  skip mRNA z-score imputation (region data incomplete)")
  } else {
    # Weighted average of available region z-scores
    num <- coalesce(df$rnafold_zscore_5utr, 0) * coalesce(df$length_5utr, 0) +
      coalesce(df$rnafold_zscore_cds,  0) * coalesce(df$length_cds,  0) +
      coalesce(df$rnafold_zscore_3utr, 0) * coalesce(df$length_3utr, 0)
    
    den <- if_else(!is.na(df$rnafold_zscore_5utr), coalesce(df$length_5utr, 0), 0) +
      if_else(!is.na(df$rnafold_zscore_cds),  coalesce(df$length_cds,  0), 0) +
      if_else(!is.na(df$rnafold_zscore_3utr), coalesce(df$length_3utr, 0), 0)
    
    df$rnafold_zscore_mrna_imputed <- if_else(den > 0, num / den, NA_real_)
    
    if ("rnafold_zscore_mrna" %in% names(df)) {
      df$rnafold_zscore_mrna <- coalesce(df$rnafold_zscore_mrna, df$rnafold_zscore_mrna_imputed)
    } else {
      df$rnafold_zscore_mrna <- df$rnafold_zscore_mrna_imputed
    }
  }
  
  if (have_sc) {
    # Score: additive across regions
    all_na <- is.na(df$rnafold_score_5utr) &
      is.na(df$rnafold_score_cds)  &
      is.na(df$rnafold_score_3utr)
    sum_score <- rowSums(cbind(df$rnafold_score_5utr,
                               df$rnafold_score_cds,
                               df$rnafold_score_3utr),
                         na.rm = TRUE)
    df$rnafold_score_mrna_imputed <- if_else(all_na, NA_real_, sum_score)
    
    if ("rnafold_score_mrna" %in% names(df)) {
      df$rnafold_score_mrna <- coalesce(df$rnafold_score_mrna, df$rnafold_score_mrna_imputed)
    } else {
      df$rnafold_score_mrna <- df$rnafold_score_mrna_imputed
    }
  }
  
  # Recompute mRNA-level expected/delta with new filled values
  # Updated: gc_mrna -> gc_content_mrna
  if (all(c("gc_content_mrna", "length_mrna", "rnafold_score_mrna") %in% names(df))) {
    df$mfe_expected_mrna <- calculate_mfe_expected(df$gc_content_mrna, df$length_mrna)
    df$mfe_delta_mrna    <- calculate_mfe_delta(df$rnafold_score_mrna, df$mfe_expected_mrna)
  }
  
  df
}

# --- MFE per nt ---
add_mfe_per_nt <- function(df) {
  regions <- c("mrna", "5utr", "cds", "3utr", "start", "stop", "last100")
  
  for (region in regions) {
    sc_col  <- paste0("rnafold_score_", region)
    len_col <- paste0("length_", region)
    out_col <- paste0("rnafold_per_nt_", region)
    
    if (!all(c(sc_col, len_col) %in% names(df))) {
      message("  skip ", out_col, " (missing ", sc_col, " or ", len_col, ")")
      next
    }
    
    len <- df[[len_col]]
    df[[out_col]] <- dplyr::if_else(len > 0, df[[sc_col]] / len, NA_real_)
  }
  df
}

# --- Junction density -------------------------------------------------------
#' Add junction-density features (junctions per kilobase)
#'
#' Reads the v3+ canonical junction-count columns `junctions_count_{region}`
#' produced by load_junctions_wide(). (Pre-fix this guarded on the old raw
#' names `n_{region}_junctions`, which the v3 loader rename retired — so the
#' guards never passed and no density column was ever created.)
add_junction_density <- function(df) {
  # 5' UTR
  if (all(c("junctions_count_5utr", "length_5utr") %in% names(df))) {
    df$junctions_density_5utr <- 1000 * df$junctions_count_5utr / df$length_5utr
  }
  # CDS
  if (all(c("junctions_count_cds", "length_cds") %in% names(df))) {
    df$junctions_density_cds <- 1000 * df$junctions_count_cds / df$length_cds
  }
  # 3' UTR
  if (all(c("junctions_count_3utr", "length_3utr") %in% names(df))) {
    df$junctions_density_3utr <- 1000 * df$junctions_count_3utr / df$length_3utr
  }

  # mRNA total: sum the regional counts. Named junctions_count_mrna (not the
  # legacy junctions_mrna) so it shares the "Junction count" metric stem with
  # its regional siblings and dodges under one dotplot tick.
  if (all(c("junctions_count_5utr", "junctions_count_cds",
            "junctions_count_3utr") %in% names(df))) {
    df$junctions_count_mrna <- coalesce(df$junctions_count_5utr, 0) +
      coalesce(df$junctions_count_cds,  0) +
      coalesce(df$junctions_count_3utr, 0)

    if ("length_mrna" %in% names(df)) {
      df$junctions_density_mrna <- 1000 * df$junctions_count_mrna / df$length_mrna
    }
  }

  df
}

#' Convert raw codon / amino-acid counts into row-normalised fractions
#'
#' Counts conflate composition with CDS length. This step divides each
#' codon count by the per-row sum of all 64 codon counts (and likewise
#' for the ~20 AA counts), so each family sums to 1.0 per transcript and
#' the columns express pure composition.
#'
#' Replaces the count columns in-place — the canonical schema now treats
#' `codon_<nnn>_cds` and `aa_<x>_cds` as fractions.
add_codon_aa_fractions <- function(df) {
  codon_cols <- grep("^codon_[acgt]{3}_cds$", names(df), value = TRUE)
  aa_cols    <- grep("^aa_[a-z]_cds$",        names(df), value = TRUE)
  
  normalise <- function(cols) {
    mat   <- as.matrix(df[, cols, drop = FALSE])
    total <- rowSums(mat, na.rm = TRUE)
    safe  <- ifelse(total > 0, total, NA_real_)
    df[, cols] <<- mat / safe
  }
  
  if (length(codon_cols) >= 2) normalise(codon_cols)
  if (length(aa_cols)    >= 2) normalise(aa_cols)
  
  df
}

# --- Drop all-NA columns ----------------------------------------------------

#' Drop columns that are entirely NA
drop_all_na_columns <- function(df) {
  keep <- !vapply(df, function(x) all(is.na(x)), logical(1))
  df[, keep, drop = FALSE]
}


# --- Orchestrator -----------------------------------------------------------

#' Apply the full feature-engineering chain to an assembled dataframe
#'
#' Order matters: junction aggregation must happen before
#' MFE expected/delta is computed so every region ends up with the full set.
#' @export
engineer_features <- function(df) {
  df |>
    add_mfe_expected_and_delta() |>
    impute_mrna_mfe() |>
    add_mfe_per_nt() |>
    add_junction_density() |>
    add_codon_aa_fractions() |>
    drop_all_na_columns()
}
