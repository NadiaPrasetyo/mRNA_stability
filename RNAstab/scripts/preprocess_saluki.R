# =============================================================================
# Preprocess Saluki predictions
# =============================================================================
# Reads the raw Saluki HDF5 prediction/target files and the accompanying
# genes.tsv, and writes a clean `saluki_predictions.rds` with two columns:
#   - gene_id
#   - saluki_prediction
#
# This is a ONE-OFF preprocessing step. The pipeline then consumes the .rds
# instead of re-parsing HDF5s on every run.
#
# Usage (from the project root):
#   Rscript scripts/preprocess_saluki.R --species human \
#       --base salukiData/deeplearning
#
# Or interactively:
#   source("scripts/preprocess_saluki.R")
#   preprocess_saluki("human", "data/raw/salukiData/deeplearning")
# =============================================================================

# Load the pipeline (for species_path, config)
source("R/load_all.R")

if (!requireNamespace("rhdf5", quietly = TRUE)) {
  stop("rhdf5 is required. Install with:\n",
       '  install.packages("BiocManager"); BiocManager::install("rhdf5")')
}


#' Load and average HDF5 predictions or targets across folds
#'
#' @param base_dir Directory containing fold_run subdirectories.
#' @param file_pattern "preds.h5" or "targets.h5".
#' @param test_set_filter Subdirectory name identifying the test set (e.g. "test0").
#' @return A dataframe with columns fold, seqnum, value.
load_and_average <- function(base_dir, file_pattern, test_set_filter) {
  all_files <- list.files(base_dir, pattern = file_pattern,
                          full.names = TRUE, recursive = TRUE)
  target_files <- all_files[grep(test_set_filter, all_files)]

  if (length(target_files) == 0) {
    stop("No files found for filter '", test_set_filter,
         "' in '", base_dir, "'.")
  }

  combined <- do.call("rbind", lapply(target_files, function(file) {
    if (grepl("preds.h5", file)) {
      vals <- t(rhdf5::h5read(file, "preds"))
    } else if (grepl("targets.h5", file)) {
      vals <- rhdf5::h5read(file, "targets")[1, 1, ]
    } else {
      return(NULL)
    }

    parts <- strsplit(file, .Platform$file.sep)[[1]]
    idx   <- which(parts == test_set_filter)
    if (length(idx) == 0 || idx == 1) return(NULL)

    fold_run <- strsplit(parts[idx - 1], "_")[[1]]
    data.frame(
      value  = as.numeric(vals),
      fold   = fold_run[1],
      run    = fold_run[2],
      seqnum = seq_along(vals)
    )
  }))

  if (is.null(combined) || nrow(combined) == 0) {
    stop("Failed to process any files.")
  }

  aggregate(value ~ fold + seqnum, data = combined, FUN = mean)
}


#' Build the gene ID map from genes.tsv, keyed by (fold, seqnum)
build_gene_id_map <- function(base_dir) {
  tsv_files <- list.files(base_dir, pattern = "genes.tsv",
                          recursive = TRUE, full.names = TRUE)
  if (length(tsv_files) == 0) {
    stop("No 'genes.tsv' found under ", base_dir)
  }

  gene_data <- read.table(tsv_files[1], header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)

  # Local seqnum that restarts per (Fold, Split)
  gene_data$seqnum <- ave(gene_data$num, gene_data$Fold, gene_data$Split,
                          FUN = seq_along)
  gene_data$fold <- paste0("f", gene_data$Fold)

  data.frame(
    fold            = gene_data$fold,
    seqnum          = gene_data$seqnum,
    gene_id         = gene_data$Gene,
    split           = gene_data$Split,
    stringsAsFactors = FALSE
  )
}


#' Run the full Saluki preprocessing
#'
#' @param species The species to write the output for (determines output path).
#' @param base_dir Directory containing the Saluki HDF5 outputs and genes.tsv.
#' @param train_dir Subdirectory of base_dir holding the actual preds/targets
#'   (default "train_gru").
#' @param test_set The test-set fold name (default "test0").
preprocess_saluki <- function(species,
                              base_dir,
                              train_dir = "train_gru",
                              test_set  = "test0") {
  data_dir <- file.path(base_dir, train_dir)

  message("Loading predictions...")
  preds_avg <- load_and_average(data_dir, "preds.h5", test_set)
  names(preds_avg)[names(preds_avg) == "value"] <- "predicted_half_life"

  message("Loading gene ID map...")
  gene_map <- build_gene_id_map(base_dir)

  merged <- merge(preds_avg, gene_map, by = c("fold", "seqnum"))

  out <- unique(merged[, c("gene_id", "predicted_half_life")])
  names(out)[2] <- "saluki_prediction"

  # Write to the species' raw data folder
  out_dir <- file.path(RAW_DIR, SPECIES_CONFIG[[species]]$dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(out_dir, SPECIES_CONFIG[[species]]$saluki_rds)

  saveRDS(out, out_path)
  message("Wrote ", nrow(out), " Saluki predictions to ", out_path)
  invisible(out_path)
}


# --- Command-line dispatch ---------------------------------------------------
# Allow running as `Rscript scripts/preprocess_saluki.R --species human --base <path>`

if (sys.nframe() == 0 && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  species <- args[which(args == "--species") + 1]
  base    <- args[which(args == "--base")    + 1]
  if (length(species) == 0 || length(base) == 0) {
    stop("Usage: Rscript scripts/preprocess_saluki.R --species <name> --base <dir>")
  }
  preprocess_saluki(species, base)
}
