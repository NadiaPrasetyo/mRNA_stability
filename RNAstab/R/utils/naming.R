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
#
# --- v4 pseudo-region retirement (CACHE_VERSION 4L) --------------------------
# v3 introduced four pseudo-region tokens — `transcript`, `window`, `core`,
# `full`. v4 removes all of them. Every region-bearing column ends in one of
# the eight real region tokens (see `REGIONS` in R/config.R).
#
# The display-formatting consequences:
#   * Whole-transcript scalars (architecture, uORF, probing) carry the real
#     `mrna` suffix instead of `transcript`. They now render with a trailing
#     "mRNA" like every other mRNA-region column — e.g.
#     `intron_length_mean_mrna` -> "Mean intron length mRNA". The old
#     `transcript` rule (which rendered to "" and was stripped) is gone.
#   * The three NMD fragility windows are collapsed to one model carrying
#     `mrna`. The NMD metric-prefix rules are unchanged; combined with the
#     `mrna` region rule they yield e.g.
#     `nmd_fragile_codon_density_mrna` -> "NMD fragile codon density mRNA".
#     The `distal`/`nearest`/`any` region rules are gone.
# =============================================================================


# Substitutions applied in order. Later rules operate on the output of earlier
# ones. Most-specific patterns first within each family. The ordering rules
# that matter:
#   * `^junctions_count_` (v3) and `^junctions_density_` MUST both precede
#     the general `^junctions_` rule (which renders as "Junctions ").
#   * NMD metric prefix rules render to "NMD <metric> "; the trailing space
#     lets the `mrna` region rule apply cleanly.
#   * Region-suffix rules are last so they see the bare region token after
#     all prefix rules have consumed the leading portion of the column name.
REPLACEMENTS <- list(
  # --- NMD Patterns (v4: single model, metric prefix, mrna region suffix) ---
  list("^nmd_fragile_codon_density_", "NMD fragile codon density "),
  list("^nmd_alt_stop_density_",      "NMD alt-stop density "),
  list("^nmd_fragile_codon_count_",   "NMD fragile codon count "),
  list("^nmd_alt_stop_count_",        "NMD alt-stop count "),  
  list("nmd_transversion_fragile_codon_density_", "NMD transversion fragile codon density "),
  list("nmd_snv_fragile_codon_density_",         "NMD fragile codon density "),
  list("nmd_alt_stop_codon_density_", "NMD alt-stop count "),
  list("nmd_transition_fraction_of_snv_fragile_", "NMD transition fragile codon fraction "),
  list("nmd_transition_fragile_codon_density_", "NMD transition fragile codon density "),

  # --- Architecture (v4: whole-transcript scalars, *_mrna suffix) ---
  list("^intron_length_mean_",         "Mean intron length "),
  list("^exon_length_first_",          "First exon length "),
  list("^exon_length_last_",           "Last exon length "),
  list("^exon_count_internal_",        "Number of internal exons "),
  list("^noncoding_length_fraction_",  "Non-coding length fraction "),

  # --- Sequence Basic ---
  list("^frac_a",                "Adenine fraction"),
  list("^frac_c",                "Cytosine fraction"),
  list("^frac_g",                "Guanine fraction"),
  list("^frac_u",                "Uracil fraction"),
  list("^at_skew_",                    "AT skew "),
  list("^gc_skew_",                    "GC skew "),
  list("^purine_ratio_",               "Purine ratio "),
  list("^amino_ratio_",                "Amino ratio "),
  list("^gc_content_",                 "GC content "),

  # --- Junctions (v3: count is feature-first; EEJ distances metric-first) ---
  list("^junctions_count_",            "Junction count "),
  list("^eej_dist_downstream",             "Closest EEJ dist. downstream"),
  list("^eej_dist_upstream_",              "Closest EEJ dist. upstream "),

  # --- uORFs (v4: whole-transcript scalars, *_mrna suffix) ---
  list("^uorf_count_",                 "uORF count "),
  list("^uorf_present_",               "Has uORF "),
  list("^dist_cap_to_first_uatg_",     "Dist. cap to first uATG "),

  # --- Identifiers / metadata ---
  list("^halflife$",                       "Half-life"),
  list("^gene_name$",                      "Gene name"),
  list("^gene_id$",                        "Ensembl gene ID"),
  list("^transcript_id$",                  "Transcript ID"),
  list("^translation_efficiency$",         "Translation efficiency"),
  list("^saluki_prediction$",              "Saluki prediction"),
  list("^prediction_difference$",          "Prediction difference"),
  list("^species$",                        "Species"),

  # --- Compound metric tokens (handled before component splitting) ---
  list("^rnafold_zscore_",                 "MFE z-score "),
  list("^rnafold_score_",                  "MFE "),
  list("^mfe_expected_",                   "MFE expected "),
  list("^mfe_delta_",                      "MFE \u0394 "),
  list("^mfe_median_",                     "MFE median "),
  list("^mfe_pval_",                       "MFE p-value "),
  list("^rnafold_per_nt_",                 "MFE/nt "),
  list("^rnalfold_zscore_",                    "Local MFE z-score "),
  list("^rnalfold_score_",                     "Local MFE "),
  list("^rnalfold_pval_",                      "Local MFE p-value "),
  list("^rnaup_zscore_",                   "RNAup z-score "),
  list("^rnaup_score_",                    "RNAup "),
  list("^rnaup_pval_",                     "RNAup p-value "),
  list("^junctions_density_",              "ORF Junction density "),
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
  list("^length_",                         "Length "),
  list("^gc_",                             "GC "),
  list("^cai$",                            "CAI"),
  list("^expression$",                     "Expression"),
  list("^orfexondensity$",                 "ORF-exon density"),
  list("^nuc_ratio_a",                     "Adenine ratio"),
  list("^nuc_ratio_c",                     "Cytosine ratio"),
  list("^nuc_ratio_g",                     "Guanine ratio"),
  list("^nuc_ratio_u",                     "Uracil ratio"),
  list("^codon_",                          "Codon freq. "),  # e.g. codon_aaa_cds -> "Codon freq. aaa CDS"
  list("^aa_",                             "AA freq. "),     # e.g. aa_l_cds      -> "AA freq. l CDS"
  list("^gini_nucleoplasm_",               "icSHAPE Gini nucleoplasm "),
  list("^gini_cytoplasm_",                 "icSHAPE Gini cytoplasm "),
  # list("^shape_",                          "icSHAPE "),
  # list("^keth_",                           "Keth-seq "),

  # --- Region suffixes ---
  # Each rule matches a single leading separator — space OR underscore —
  # plus the region token, end-anchored (`[ _]<region>$`). The separator is
  # reproduced as a space in the replacement.
  #
  # Why match the separator explicitly rather than `\b<region>$`: the region
  # is always preceded by `_` in a raw column, and `_` is a regex *word*
  # character, so `\b` sits at no boundary there. `\bmrna$` only fired once
  # an earlier prefix rule had already turned that `_` into a space — it
  # silently failed on columns where a prefix rule left a residual `_`
  # before the region (e.g. `shape_score_mrna`, `exon_length_first_fraction_
  # mrna`), leaking a lowercase token into the display string.
  #
  # End-anchoring plus a literal `[ _]` separator also correctly avoids the
  # `alt-stop` collision: a hyphen is neither space nor underscore, so
  # `[ _]stop$` cannot fire on the `stop` inside `alt-stop` (whereas the old
  # `\bstop$` would have — a hyphen is a non-word char, so a boundary sits
  # before `stop`).
  list("[ _]5utr$",                        " 5' UTR"),
  list("[ _]3utr$",                        " 3' UTR"),
  list("[ _]cds$",                         " CDS"),
  list("[ _]mrna$",                        " mRNA"),
  list("[ _]utrpair$",                     " UTR interactions"),
  list("[ _]last100$",                     " last 100 nt"),
  list("[ _]start$",                       " start codon"),
  list("[ _]stop$",                        " stop codon")
)


#' Format a canonical column name into a display string
#'
#' @param col_name Character (length-1 or vector). A canonical column name.
#' @return Character vector of the same length.
#' @examples
#' format_col_name("rnafold_zscore_5utr")   # "MFE z-score 5' UTR"
#' format_col_name("length_cds")            # "Length CDS"
#' format_col_name("halflife")              # "Half-life"
#' # v4 examples:
#' format_col_name("junctions_count_5utr")            # "Junction count 5' UTR"
#' format_col_name("eej_dist_first_start")            # "First EEJ dist. from start codon"
#' format_col_name("intron_length_mean_mrna")         # "Mean intron length mRNA"
#' format_col_name("nmd_fragile_codon_density_mrna")  # "NMD fragile codon density mRNA"
#' format_col_name("uorf_count_mrna")                 # "uORF count mRNA"
#' @export
format_col_name <- function(col_name) {
  vapply(col_name, format_single_name, character(1), USE.NAMES = FALSE)
}

format_single_name <- function(name) {
  # Codon / amino-acid composition: token uppercased, region suffix dropped.
  # Cheap special-case avoids 84 exact-match rules in REPLACEMENTS and
  # avoids needing PCRE2 case-folding (which R's sub() doesn't support).
  m <- regmatches(name, regexec("^codon_([acgtu]{3})_cds$", name))[[1]]
  if (length(m) == 2) return(paste0("Codon freq. ", toupper(m[2])))
  
  m <- regmatches(name, regexec("^aa_([a-z])_cds$", name))[[1]]
  if (length(m) == 2) return(paste0("AA freq. ", toupper(m[2])))
  
  for (rule in REPLACEMENTS) {
    name <- sub(rule[[1]], rule[[2]], name)
  }
  name <- gsub("_", " ", name, fixed = TRUE)
  name <- gsub("  +", " ", name)
  trimws(name)
}
