# =============================================================================
# Raw data loaders
# =============================================================================
# One function per source file. Every loader:
#   * takes a species argument
#   * returns a tibble with canonical column names (lowercase, snake_case)
#   * returns NULL if the source file is missing (pipeline tolerates this)
#   * normalises the `region` vocabulary to lowercase
#   * ensures the region token (real or pseudo) is the LAST token of any
#     region-bearing column — see the v3 migration notes below.
#
# Region-level loaders return LONG-form data (one row per transcript × region).
# The pipeline pivots to wide form in `pipeline/assemble.R`.
#
#   load_junctions_wide : n_<region>_junctions   -> junctions_count_<region>
#                         <anchor>_dist_<pos>_eej -> eej_dist_<pos>_<anchor>
#   load_architecture   : intron_mean / *_exon_* -> *_mrna
#   load_uorfs          : n_uorfs / has_uorf ...  -> uorf_*_mrna
#   load_nmd_fragility  : single model, *_metric -> nmd_<metric>_mrna
#
# If a loader here doesn't match the shape of your raw file, update the loader
# rather than the pipeline — loaders are the stable boundary.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})


# --- Helpers -----------------------------------------------------------------

normalise_region <- function(df) {
  if ("region" %in% names(df)) {
    df$region <- tolower(df$region)
  } else if ("Sequence_Type" %in% names(df)) {
    df <- dplyr::rename(df, region = Sequence_Type)
    df$region <- tolower(df$region)
  }
  
  if ("region" %in% names(df)) {
    # Apply alias table.
    hits <- df$region %in% names(REGION_ALIASES)
    if (any(hits)) df$region[hits] <- unname(REGION_ALIASES[df$region[hits]])
    
    # Validation safety net: anything unrecognised after aliasing is a bug
    # waiting to happen — it'll silently fail the region-suffix-last
    # invariant once it gets pivoted to wide form.
    unknown <- setdiff(unique(df$region), REGIONS)
    if (length(unknown) > 0) {
      warning("normalise_region: unknown region token(s): ",
              paste(unknown, collapse = ", "),
              " (expected one of: ", paste(REGIONS, collapse = ", "), ")")
    }
  }
  
  df
}


# Rename the gene_id column
rename_gene_id <- function(df) {
  if ("ensembl_gene_id" %in% names(df)) {
    df <- dplyr::rename(df, gene_id = ensembl_gene_id)
  }
  df
}

# Rename the transcript column
rename_transcript <- function(df) {
  if ("tx_id" %in% names(df)) {
    df <- dplyr::rename(df, transcript_id = tx_id)
  }
  df
}

# Standardise column names to lowercase. Useful after reading files with
# inconsistent casing.
lowercase_names <- function(df) {
  names(df) <- tolower(names(df))
  df
}

# Rename columns according to a mapping, skipping any not present in df.
# `mapping` is a named character vector in dplyr::rename form: c(new = "old").
# Columns named in the mapping but absent from df are silently ignored, which
# keeps loaders tolerant of partial / incomplete raw files.
rename_if_present <- function(df, mapping) {
  present <- mapping[mapping %in% names(df)]
  if (length(present) > 0) df <- dplyr::rename(df, !!!present)
  df
}


# Read a CSV and return NULL if it doesn't exist (with a quiet message).
# Any file-level failure surfaces as NULL so pipeline can carry on.
# Updated to handle both .csv and .tsv based on extension
read_if_exists <- function(path, ...) {
  if (!file.exists(path)) {
    message("  skip (missing): ", path)
    return(NULL)
  }

  ext <- tools::file_ext(path)
  if (ext == "tsv") {
    readr::read_tsv(path, show_col_types = FALSE, ...) |> suppressWarnings()
  } else {
    readr::read_csv(path, show_col_types = FALSE, ...) |> suppressWarnings()
  }
}

# --- Metadata loaders --------------------------------------------------------

#' Load transcript-to-gene mapping.
#' Expected columns: transcript_id, gene_id, gene_name.
load_transcripts <- function(species) {
  df <- read_if_exists(species_path(species, "transcripts.csv"))
  if (is.null(df)) return(NULL)
  df |>
    lowercase_names() |>
    dplyr::select(dplyr::any_of(c("transcript_id", "gene_id", "gene_name")))
}


#' Load half-life data.
#' Expected columns: gene_id (or gene_name), halflife.
load_halflife <- function(species) {
  df <- read_if_exists(species_path(species, "halflife.csv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names()
}


#' Load translation efficiency data.
#' Expected columns: gene_name (or gene_id), translation_efficiency.
load_translation_efficiency <- function(species) {
  df <- read_if_exists(species_path(species, "translation_efficiency.csv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names() |>
    tidyr::separate_wider_delim(
      cols = gene_id, delim = ".",
      names = c("gene_id", "gene_id_version"),
      too_many = "merge", too_few = "align_start"
    ) |> dplyr:: rename(translation_efficiency = mean_te) |>
    dplyr::select(gene_id, translation_efficiency)
}


#' Load pre-processed Saluki predictions (see scripts/preprocess_saluki.R).
#' Expected columns: gene_id, saluki_prediction.
load_saluki_predictions <- function(species) {
  cfg <- SPECIES_CONFIG[[species]]
  if (is.null(cfg$saluki_rds)) return(NULL)
  path <- species_path(species, cfg$saluki_rds)
  if (!file.exists(path)) {
    message("  skip (missing): ", path)
    return(NULL)
  }
  readRDS(path)
}


# --- Regional (long-form) loaders -------------------------------------------
# These all return long-form data keyed by (transcript_id, region).

#' Load RNAfold data (global secondary structure).
#' Columns: transcript_id, region, score, zscore, median, pval.
load_rnafold <- function(species) {
  df <- read_if_exists(species_path(species, "rnafold_results.csv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names() |> normalise_region()
}


#' Load RNALfold data (local secondary structure).
#' Columns: transcript_id, region, rnalfold_score, rnalfold_zscore, rnalfold_pval.
load_rnalfold <- function(species) {
  df <- read_if_exists(species_path(species, "rnalfold_results.csv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names() |> normalise_region()
}


#' Load RNAup data (RNA-RNA interaction).
#' Columns: transcript_id, region, rnaup_score, rnaup_zscore, rnaup_pval.
load_rnaup <- function(species) {
  df <- read_if_exists(species_path(species, "rnaup_results.csv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names() |> normalise_region()
}


#' Load sequence basic metrics (Long Form)
#' Replaces old load_length_gc and load_skews.
load_sequence_basic <- function(species) {
  df <- read_if_exists(species_path(species, "sequence_basic.tsv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names() |> normalise_region() |>
    dplyr::select(-dplyr::any_of("gene_id")) # gene_id handled by transcripts
}

#' Load NMD fragility model (Wide Form)
load_nmd_fragility <- function(species) {
  df <- read_if_exists(species_path(species, "nmd_fragility.tsv"))
  if (is.null(df)) return(NULL)

  df |>
    lowercase_names() |>
    dplyr::select(-dplyr::any_of(c("gene_id", "strand", "model", "cds_length", "zone_length",
                                   "n_transition_fragile_codons", "n_transversion_fragile_codons",
                                   "n_snv_fragile_codons", "n_alt_stop_codons"
                                   ))) |>
    # Normalise count metric stems to "feature first" form.
    rename_if_present(c(
      fragile_codon_count      = "n_fragile_codons",
      alt_stop_count           = "n_alt_stops"
    )) |>
    # Prefix with nmd_, suffix with the mrna region token.
    affix_payload(prefix = "nmd_", suffix = "_mrna", keys = "transcript_id")
}

#' Load Architecture features (Wide Form)
#'
#' v4: architecture metrics are whole-transcript scalars (one value per
#' transcript). They carry the real `mrna` region suffix and "feature first"
#' stems so they obey the region-suffix-last invariant and group cleanly
#' under FEATURE_PATTERNS. Note `intron_length_mean_mrna` is a mild abuse —
#' introns are spliced out of the mRNA — accepted in v4 to keep one region
#' vocabulary; the metric stem still reads correctly.
load_architecture <- function(species) {
  df <- read_if_exists(species_path(species, "architecture.tsv"))
  if (is.null(df)) return(NULL)
  df |>
    lowercase_names() |>
    dplyr::select(-dplyr::any_of(c("gene_id", "strand", "n_exons", "first_exon_length"))) |>
    rename_if_present(c(
      intron_length_mean_mrna        = "intron_mean",
      exon_length_first_mrna         = "first_exon_length",
      exon_length_last_mrna          = "last_exon_length",
      exon_count_internal_mrna       = "n_internal_exons",
      noncoding_length_fraction_mrna = "length_noncoding_fraction"
    ))
}

#' Load uORF features (Updated for new uorf.tsv schema)
#'
#' v4: the whole-transcript uORF scalars (count, presence flag, cap-to-uATG
#' distance) carry the real `mrna` region suffix. Any other `uorf_*` columns
#' in the raw file pass through unchanged — review them if they should also
#' carry a region token.
load_uorfs <- function(species) {
  df <- read_if_exists(species_path(species, "uorf.tsv"))
  if (is.null(df)) return(NULL)
  df |>
    lowercase_names() |>
    dplyr::select(-dplyr::any_of("gene_id")) |>
    rename_if_present(c(
      uorf_count_mrna             = "n_uorfs",
      uorf_present_mrna           = "has_uorf",
      dist_cap_to_first_uatg_mrna = "dist_cap_to_first_uatg"
    ))
}

#' Load Codon and Amino Acid counts (Wide Form)
#'
#' Raw columns are bare codon triplets / single-letter AA codes
#' (`codon_aaa`, `aa_l`, …) and are CDS-only. Per the region-suffix-last
#' invariant the loader stamps `_cds` on every payload column. The counts
#' are converted to row-normalised fractions in add_codon_aa_fractions().
load_codon_aa_counts <- function(species) {
  df <- read_if_exists(species_path(species, "codon_aa_counts.tsv"))
  if (is.null(df)) return(NULL)
  df |>
    lowercase_names() |>
    dplyr::select(-dplyr::any_of(c("gene_id", "strand"))) |>
    affix_payload(suffix = "_cds", keys = "transcript_id")
}

#' Load Junctions and Distances (Wide Form)
load_junctions_wide <- function(species) {
  df <- read_if_exists(species_path(species, "junctions.tsv"))
  if (is.null(df)) return(NULL)
  df |>
    lowercase_names() |>
    dplyr::select(-dplyr::any_of(c("gene_id", "strand"))) |>
    rename_if_present(c(
      # Junction counts: feature first, region last.
      junctions_count_5utr = "n_5utr_junctions",
      junctions_count_cds  = "n_cds_junctions",
      junctions_count_3utr = "n_3utr_junctions",
      # EEJ distances: metric first, anchor (start/stop) region last.
      eej_dist_upstream_start = "start_dist_closest_upstream",
      eej_dist_downstream_start  = "start_dist_closest_upstream",
      eej_dist_upstream_stop  = "stop_dist_closest_upstream",
      eej_dist_downstream_stop   = "stop_dist_closest_downstream"
    ))
}

#' Load stop-free region sizes.
#' Columns: transcript_id, region, stopfree.
load_stopfree <- function(species) {
  df <- read_if_exists(species_path(species, "stopfree.tsv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names() |> normalise_region()
}

# --- Wide-form transcript-level loaders --------------------------------------
# These return one row per transcript, no region column.

#' Load Codon Adaptation Index.
#' Expected columns: transcript_id, cai.
load_cai <- function(species) {
  df <- read_if_exists(species_path(species, "cai.csv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names()
}


#' Load nucleotide translation efficiency features (whatever your upstream
#' `nte` step produces). Expected columns: transcript_id plus nte_* columns.
load_nte <- function(species) {
  df <- read_if_exists(species_path(species, "nte.csv"))
  if (is.null(df)) return(NULL)
  df |> lowercase_names()
}


#' Load Agarwal features (shared file, subset by gene_id).
#' Expected columns: gene_id (or ENSID), expression, orfexondensity.
load_agarwal_features <- function(species) {
  df <- read_if_exists(species_path(species, "agarwal_features.csv"))
  if (is.null(df)) return(NULL)
  df <- lowercase_names(df)
  if ("ensid" %in% names(df)) df <- dplyr::rename(df, gene_id = ensid)
  df |> dplyr::select(dplyr::any_of(c("gene_id", "orfexondensity")))
}

#' Load Gini probing scores (icSHAPE-derived) from a shared long-form file.
#'
#' Replaces the former load_dani_probing. The shared file `gini_long.tsv` holds
#' both species in a `species` column and is LONG over two dimensions:
#' `compartment` (nucleoplasm / cytoplasm) and `region` (mrna / 5utr / cds /
#' 3utr), one `gini` value per (gene_id, compartment, region).
#'
#' Keyed on `gene_id`, NOT transcript_id: the icSHAPE data came from an older
#' annotation whose transcript IDs don't match the rest of the pipeline, but
#' the gene_ids do (and are unversioned, matching transcripts.csv). The new
#' gini pipeline collapsed to one transcript per gene, so each
#' (gene_id, compartment, region) cell is a single value — the pivot below is
#' guaranteed not to produce list-columns. `transcript_id` and `n_valid` are
#' dropped.
#'
#' Returns wide form keyed on gene_id, with columns `gini_<compartment>_<region>`
#' (region token last, per the schema invariant). Routed through the gene-level
#' join in build_dataset(); gini is therefore broadcast to every transcript of
#' a gene.
load_gini_probing <- function(species) {
  df <- read_if_exists(shared_path("gini_long.tsv"))
  if (is.null(df)) return(NULL)
  
  df <- df |>
    lowercase_names() |>      # MRNA -> mrna handled by normalise_region below
    normalise_region()        # lowercases the region column
  
  if (!"species" %in% names(df)) return(NULL)
  df <- dplyr::filter(df, species == !!species)
  if (nrow(df) == 0) return(NULL)
  
  df |>
    dplyr::select(dplyr::any_of(c("gene_id", "compartment", "region", "gini"))) |>
    tidyr::pivot_wider(
      id_cols     = gene_id,
      names_from  = c(compartment, region),
      values_from = gini,
      names_glue  = "gini_{compartment}_{region}"
    )
}

#' #' Load Dani probing data (icSHAPE + Keth-seq) from a shared file.
#' #'
#' #' The shared file holds data for multiple genomes in genome-suffixed columns
#' #' (e.g. `icSHAPE_hg38_*` vs `icSHAPE_mm10_*`). This picks the right subset
#' #' for the requested species using SPECIES_CONFIG[[species]]$dani patterns,
#' #' and renames them to generic `shape_*` / `keth_*` columns.
#' load_dani_probing <- function(species) {
#'   df <- read_if_exists(shared_path("dani_probing.csv"))
#'   if (is.null(df)) return(NULL)
#' 
#'   cfg <- SPECIES_CONFIG[[species]]$dani
#'   id_cols <- c("gene_name", "transcript_id", "gene_id")
#'   id_cols <- intersect(id_cols, names(df))
#' 
#'   # Pull id columns + shape/keth columns for this species + any extras
#'   keep <- c(
#'     id_cols,
#'     grep(cfg$shape_pattern, names(df), value = TRUE),
#'     grep(cfg$keth_pattern,  names(df), value = TRUE),
#'     intersect(cfg$extra_cols, names(df))
#'   )
#'   df <- df[, unique(keep), drop = FALSE]
#' 
#'   # Rename: icSHAPE_hg38_* → shape_*, Keth_seq_hg38_* → keth_*
#'   names(df) <- stringr::str_replace(names(df), cfg$shape_pattern, "shape")
#'   names(df) <- stringr::str_replace(names(df), cfg$keth_pattern,  "keth")
#'   names(df) <- tolower(names(df))
#'   # Drop any columns containing the string "cvg"
#'   df <- df[, !grepl("cvg", names(df)), drop = FALSE]
#'   df <- affix_payload(df, suffix = "_mrna", keys = c("gene_name", "pairwise_identy_CDS"))
#'   df
#' }


# --- Internal helpers --------------------------------------------------------

# Reduce a list of optional (possibly NULL) dataframes by full_join on an id.
Reduce_safely_by_id <- function(dfs, id_col) {
  dfs <- Filter(Negate(is.null), dfs)
  if (length(dfs) == 0) return(NULL)
  Reduce(function(a, b) dplyr::full_join(a, b, by = id_col), dfs)
}
