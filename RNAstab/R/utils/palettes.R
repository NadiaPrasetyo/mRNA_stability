# =============================================================================
# Feature group palettes and shape mappings
# =============================================================================
# Hand-curated colour scheme for FEATURE_PATTERNS keys, organised by the
# SUPERGROUPS membership defined in R/config.R:
#
#   structure  — purples, magentas, oranges, yellow
#   intrinsic  — blues, greens, teals
#   splicing   — reds, browns
#   regulatory — distinct purples (avoids structure family clash)
#   other      — greys (standalone columns and catch-all)
#
# Plus Okabe-Ito colourblind-friendly palette for REGIONS, used when region
# is the visual variable rather than the grouping (see feature_correlation_
# dotplot vs feature_response_scatter for the two patterns).
#
# v4: the v3 pseudo-region tokens (`transcript`, `window`, `core`, `full`)
# have been retired. Whole-transcript scalars and the single NMD model now
# use the real `mrna` region, so REGION_COLOURS/SHAPES/DISPLAYS cover only
# the eight real regions plus the `none` fallback.
#
# Adapted from the legacy plotting scripts. To override: pass `palette = ...`
# or `region_colours = ...` to plot functions. To extend: add keys here.
# =============================================================================


# -----------------------------------------------------------------------------
# Feature-group colours (used when GROUP is the visual variable)
# -----------------------------------------------------------------------------

#' Colour for each FEATURE_PATTERNS group plus a few standalone columns.
#' @export
FEATURE_GROUP_COLOURS <- c(

  # --- Structure: purples + oranges + yellow ---
  rnafold_scores         = "#2D004B",  # very dark purple
  rnafold_zscores        = "#E7298A",  # vivid magenta
  rnafold_per_nt         = "#D4115A",
  mfe_deltas             = "#542788",  # deep violet
  mfe_expected           = "#8073AC",  # muted purple
  rnalfold_scores            = "#FF7F00",  # vivid orange
  rnalfold_zscores           = "#B35806",  # burnt orange
  rnaup                  = "#FDB863",  # light orange
  probing                = "#FFD92F",  # bright yellow

  # --- Intrinsic: blues + greens + teals ---
  lengths                = "#002642",  # midnight blue
  gc                     = "#1F78B4",  # vivid blue
  stopfree               = "#A6CEE3",  # pale blue
  codon_freqs            = "#33A02C",  # vivid green
  aa_freqs               = "#B2DF8A",  # pale green
  nuc_ratios             = "#01665E",  # dark teal
  skews                  = "#80CDC1",  # light teal

  # --- Splicing: reds + browns ---
  junctions              = "#E31A1C",  # vivid red
  distances              = "#FB9A99",  # pinkish red
  architecture           = "#8C510A",  # earth brown
  nmd                    = "#B15928",  # rust

  # --- Regulatory: distinct purples (separate family from structure) ---
  uorfs                  = "#6A3D9A",  # deep purple
  orfs                   = "#CAB2D6",  # light purple

  # --- Standalone columns: greys and black ---
  cai                    = "#000000",  # pure black
  translation_efficiency = "#525252",  # dark grey
  expression             = "#7F7F7F",  # medium grey
  orfexondensity         = "#9E9E9E",  # mid-light grey
  other                  = "#C7C7C7"   # light grey (catch-all)
)


# -----------------------------------------------------------------------------
# Region colours and shapes (used when REGION is the visual variable)
# -----------------------------------------------------------------------------
# Okabe-Ito colourblind-friendly palette for the eight real regions. Order
# intentionally follows the 5' -> 3' traversal of the transcript so the
# legend reads naturally and intersect()-based ordering preserves it.

#' Colours for each REGIONS token plus a fallback.
#' @export
REGION_COLOURS <- c(
  `5utr`     = "#009E73",  # bluish green
  start      = "#006442",  # green
  cds        = "#56B4E9",  # sky blue
  mrna       = "#0072B2",  # blue
  stop       = "#E69F00",  # orange
  `3utr`     = "#D55E00",  # vermillion
  last100    = "#CC79A7",  # reddish purple
  utrpair    = "#000000",  # black
  none       = "#999999"   # grey (unrecognised / no region)
)


#' Shape codes for each REGIONS token plus a fallback.
#' @export
REGION_SHAPES <- c(
  `5utr`     = 17,   # filled triangle (up)
  start      = 2,    # open triangle (up)
  cds        = 15,   # filled square
  mrna       = 16,   # filled circle
  stop       = 6,    # open triangle (down)
  `3utr`     = 18,   # filled diamond
  last100    = 8,    # asterisk
  utrpair    = 11,   # crossed square
  none       = 4     # x (unrecognised / no region)
)


#' Display strings for each REGIONS token.
#' MUST stay in sync with the Region tokens table in PIPELINE_GUIDE §2.4 and
#' with the region rules in R/utils/naming.R. Used by format_metric_name()
#' to strip the region from a formatted column name.
#' @export
REGION_DISPLAYS <- c(
  `5utr`     = "5' UTR",
  cds        = "CDS",
  `3utr`     = "3' UTR",
  mrna       = "mRNA",
  utrpair    = "UTR interactions",
  last100    = "last 100 nt",
  start      = "start codon",
  stop       = "stop codon",
  none       = ""
)


# -----------------------------------------------------------------------------
# Display formatters
# -----------------------------------------------------------------------------

FEATURE_GROUP_DISPLAY_NAMES <- c(
  rnafold_scores         = "MFE",
  rnafold_zscores        = "MFE z-score",
  rnafold_per_nt         = "MFE/nt",  
  mfe_deltas             = "MFE Delta",
  mfe_expected           = "MFE expected",
  rnalfold_scores            = "Local MFE",
  rnalfold_zscores           = "Local MFE z-score",
  rnaup                  = "RNAup",
  probing                = "Probing",
  lengths                = "Length",
  gc                     = "GC content",
  stopfree               = "Stop-free",
  codon_freqs            = "Codon frequency",
  aa_freqs               = "Amino acid frequency",
  nuc_ratios             = "Nucleotide ratios",
  skews                  = "Skews",
  junctions              = "Junctions",
  distances              = "EEJ distance",
  architecture           = "Architecture",
  nmd                    = "NMD fragility",
  uorfs                  = "uORFs",
  orfs                   = "ORFs",
  cai                    = "CAI",
  translation_efficiency = "Translation efficiency",
  expression             = "Expression",
  orfexondensity         = "ORF-exon density",
  other                  = "Other"
)


#' Display name for a feature group key.
#'
#' Falls back to title-cased underscore replacement for keys not in the
#' lookup table. Vectorised.
#' @export
format_group_name <- function(g) {
  out <- FEATURE_GROUP_DISPLAY_NAMES[g]
  missing_idx <- is.na(out)
  if (any(missing_idx)) {
    out[missing_idx] <- tools::toTitleCase(gsub("_", " ", g[missing_idx]))
  }
  unname(out)
}


#' Display name for a column with its region suffix stripped.
#'
#' Useful when region is encoded as a separate visual (colour, shape) and
#' shouldn't be repeated in the axis label. Strategy: format the full column
#' via format_col_name(), then strip the formatted region from the end.
#' Falls back to format_col_name() unchanged for columns without a region
#' suffix.
#'
#' @param col Character vector of column names.
#' @return Character vector of display strings.
#' @examples
#' format_metric_name("length_cds")                    # "Length"
#' format_metric_name("rnafold_zscore_5utr")            # "MFE z-score"
#' format_metric_name("cai")                            # "CAI"
#' format_metric_name("intron_length_mean_mrna")        # "Mean intron length"
#' format_metric_name("nmd_fragile_codon_count_mrna")   # "NMD fragile codon count"
#' @export
format_metric_name <- function(col) {
  vapply(col, function(co) {
    tokens <- strsplit(co, "_", fixed = TRUE)[[1]]
    if (length(tokens) <= 1) return(format_col_name(co))

    last <- tokens[length(tokens)]
    if (!last %in% REGIONS) return(format_col_name(co))

    full   <- format_col_name(co)
    suffix <- paste0(" ", REGION_DISPLAYS[[last]])

    # Use endsWith + substr to avoid regex escaping. If the formatted name
    # doesn't actually end in the region display string, return it unchanged.
    if (nchar(suffix) > 1 && endsWith(full, suffix)) {
      trimws(substr(full, 1, nchar(full) - nchar(suffix)))
    } else {
      full
    }
  }, character(1), USE.NAMES = FALSE)
}


# -----------------------------------------------------------------------------
# Convenience accessors
# -----------------------------------------------------------------------------

#' Get colours for a vector of group keys, falling back to "other" grey.
#' @export
feature_colour <- function(group) {
  out <- FEATURE_GROUP_COLOURS[group]
  out[is.na(out)] <- FEATURE_GROUP_COLOURS[["other"]]
  unname(out)
}

#' Get colours for a vector of region tokens, falling back to "none" grey.
#' @export
region_colour <- function(region) {
  region <- ifelse(is.na(region) | region == "", "none", region)
  out <- REGION_COLOURS[region]
  out[is.na(out)] <- REGION_COLOURS[["none"]]
  unname(out)
}

#' Get shapes for a vector of region tokens, falling back to "none".
#' @export
region_shape <- function(region) {
  region <- ifelse(is.na(region) | region == "", "none", region)
  out <- REGION_SHAPES[region]
  out[is.na(out)] <- REGION_SHAPES[["none"]]
  unname(out)
}
