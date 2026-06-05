# =============================================================================
# Build the mouse dataset
# =============================================================================
# Same as build_human.R but for mouse data. This assumes the mouse raw files
# have been placed under data/raw/mouse/ following the same layout as human.
#
# Run:
#   Rscript scripts/build_mouse.R
#   Rscript scripts/build_mouse.R --rebuild
# =============================================================================

source("R/load_all.R")

args <- commandArgs(trailingOnly = TRUE)
rebuild <- "--rebuild" %in% args

df <- build_dataset("mouse", rebuild = rebuild)

cat("\nBuilt mouse dataset:\n")
cat("  rows:    ", nrow(df), "\n")
cat("  columns: ", ncol(df), "\n")
cat("  complete cases: ", sum(complete.cases(df)), "\n")
