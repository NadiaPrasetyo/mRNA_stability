# =============================================================================
# Build the human dataset
# =============================================================================
# Run this once to populate the cache:
#   Rscript scripts/build_human.R            # uses cache if present
#   Rscript scripts/build_human.R --rebuild  # force rebuild
# =============================================================================

source("R/load_all.R")

args <- commandArgs(trailingOnly = TRUE)
rebuild <- "--rebuild" %in% args

df <- build_dataset("human", rebuild = rebuild)

cat("\nBuilt human dataset:\n")
cat("  rows:    ", nrow(df), "\n")
cat("  columns: ", ncol(df), "\n")
cat("  complete cases: ", sum(complete.cases(df)), "\n")
