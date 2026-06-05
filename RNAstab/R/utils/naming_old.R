# =============================================================================
# Display-name formatter for plot labels
# =============================================================================
# Converts canonical column names (lowercase, underscore-separated) into
# human-readable strings suitable for plot axes, legends and titles.
#
# The pipeline stores data with a canonical schema — e.g. "rnafold_zscore_5utr",
# "length_cds". This function is the bridge between that schema and the
# pretty labels end-users see on figures.
#
# Add to `REPLACEMENTS` to teach it new patterns.
# =============================================================================


# Substitutions applied in order. Later rules operate on the output of earlier
# ones. Most-specific patterns first.
REPLACEMENTS <- list(
  # --- NMD Patterns ---
  list("^nmd_window_",             "NMD (window) "),
  list("^nmd_core_",               "NMD (core) "),
  list("^nmd_all_",                "NMD (all) "),
  list("_fragile_codon_density$",  " fragile codon density"),
  list("_alt_stop_density$",       " alt-stop density"),
  list("_n_fragile_codons$",       " fragile codon count"),
  list("_n_alt_stops$",            " alt-stop count"),
  
  # --- Architecture ---
  list("^n_internal_exons$",       "Number of internal exons"),
  list("^first_exon_length$",      "First exon length"),
  list("^last_exon_length$",       "Last exon length"),
  list("^intron_mean$",            "Mean intron length"),
  
  # --- Sequence Basic ---
  list("^frac_",                   "Fraction "),
  list("^at_skew_",                "AT skew "),
  list("^gc_skew_",                "GC skew "),
  list("^purine_ratio_",           "Purine ratio "),
  list("^amino_ratio_",            "Amino ratio "),
  list("^gc_content_",             "GC content "),
  
  # --- Junctions (New) ---
  list("^n_5utr_junctions$",       "5' UTR junction count"),
  list("^n_cds_junctions$",        "CDS junction count"),
  list("^n_3utr_junctions$",       "3' UTR junction count"),
  list("^stop_dist_",              "Dist. from stop to "),
  list("^start_dist_",             "Dist. from start to "),
  
  # --- uORFs (New) ---
  list("^n_uorfs$",                "uORF count"),
  list("^has_uorf$",               "Has uORF"),
  list("^dist_cap_to_first_uatg$", "Dist. cap to first uATG"),
  
  # --- Identifiers / metadata ---
  list("^halflife$",                       "Half-life"),
  list("^gene_name$",                      "Gene name"),
  list("^gene_id$",                "Ensembl gene ID"),
  list("^transcript_id$",                  "Transcript ID"),
  list("^translation_efficiency$",         "Translation efficiency"),
  list("^saluki_prediction$",              "Saluki prediction"),
  list("^prediction_difference$",          "Prediction difference"),
  list("^species$",                        "Species"),

  # --- Compound metric tokens (handled before component splitting) ---
  list("^rnafold_zscore_",                     "MFE z-score "),
  list("^rnafold_score_",                      "MFE "),
  list("^mfe_expected_",                   "MFE expected "),
  list("^mfe_delta_",                      "MFE \u0394 "),
  list("^mfe_median_",                     "MFE median "),
  list("^mfe_pval_",                       "MFE p-value "),
  list("^rnalfold_zscore_",                    "Local MFE z-score "),
  list("^rnalfold_score_",                     "Local MFE "),
  list("^rnalfold_pval_",                      "Local MFE p-value "),
  list("^rnaup_zscore_",                   "RNAup z-score "),
  list("^rnaup_score_",                    "RNAup "),
  list("^rnaup_pval_",                     "RNAup p-value "),
  list("^junctions_density_",              "Junction density "),
  list("^junctions_",                      "Junctions "),
  list("^orf_percent_length_",             "uORF % length "),
  list("^orf_number_",                     "uORF count "),
  list("^orf_length_",                     "uORF length "),
  list("^orfj_density$",                   "ORF-J density"),
  list("^stopfree_",                       "Stop-free "),
  list("^distance_junction_up_start$",     "Junction upstream distance (start)"),
  list("^distance_junction_down_start$",   "Junction downstream distance (start)"),
  list("^distance_junction_up_stop$",      "Junction upstream distance (stop)"),
  list("^distance_junction_down_stop$",    "Junction downstream distance (stop)"),
  list("^length_noncoding_fraction$",      "Non-coding length fraction"),
  list("^length_",                         "Length "),
  list("^gc_",                             "GC "),
  list("^cai$",                            "CAI"),
  list("^expression$",                     "Expression"),
  list("^orfexondensity$",                 "ORF-exon density"),
  list("^nuc_ratio_a",                     "Adenine ratio"),
  list("^nuc_ratio_c",                     "Cytosine ratio"),
  list("^nuc_ratio_g",                     "Guanine ratio"),
  list("^nuc_ratio_u",                     "Uracil ratio"),
  list("^freq_cds_",                       "Codon freq. CDS "),
  list("^freq_last10_",                    "Codon freq. last-10 "),
  list("^aa_freq_",                        "AA freq. "),
  list("^shape_",                          "icSHAPE "),
  list("^keth_",                           "Keth-seq "),

  # --- Region suffixes ---
  list("5utr\\b",                          "5' UTR"),
  list("3utr\\b",                          "3' UTR"),
  list("\\bcds\\b",                        "CDS"),
  list("\\bmrna\\b",                       "mRNA"),
  list("\\butrpair\\b",                      "UTR interactions"),
  list("\\blast100\\b",                    "last 100 nt"),
  list("\\bstart\\b",                      "start codon"),
  list("\\bstop\\b",                       "stop codon")
)


#' Format a canonical column name into a display string
#'
#' @param col_name Character (length-1 or vector). A canonical column name.
#' @return Character vector of the same length.
#' @examples
#' format_col_name("rnafold_zscore_5utr")   # "MFE z-score 5' UTR"
#' format_col_name("length_cds")        # "Length CDS"
#' format_col_name("halflife")          # "Half-life"
#' @export
format_col_name <- function(col_name) {
  vapply(col_name, format_single_name, character(1), USE.NAMES = FALSE)
}

format_single_name <- function(name) {
  for (rule in REPLACEMENTS) {
    name <- sub(rule[[1]], rule[[2]], name)
  }
  # Final tidy-up: convert any leftover underscores, trim whitespace
  name <- gsub("_", " ", name, fixed = TRUE)
  name <- trimws(name)
  name
}
