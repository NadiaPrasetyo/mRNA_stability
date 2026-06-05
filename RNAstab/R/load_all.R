# =============================================================================
# load_all.R
# =============================================================================
# Source this file to load the entire pipeline into an R session:
#
#   source("R/load_all.R")
#   df <- build_dataset("human")
#
# Files are sourced in dependency order:
#   1. config.R         (constants, paths, species registry)
#   2. utils/           (no dependencies on other pipeline files)
#   3. io/              (depends on config)
#   4. features/        (depends on utils + config)
#   5. pipeline/        (depends on everything)
# =============================================================================

.PIPELINE_ROOT <- if (exists(".PIPELINE_ROOT")) .PIPELINE_ROOT else "R"

.source_dir <- function(subdir) {
  files <- list.files(file.path(.PIPELINE_ROOT, subdir),
                      pattern = "\\.R$", full.names = TRUE)
  for (f in files) source(f, local = FALSE)
}

source(file.path(.PIPELINE_ROOT, "config.R"))
.source_dir("utils")
.source_dir("io")
.source_dir("features")
.source_dir("pipeline")

message("Pipeline loaded. Available species: ",
        paste(names(SPECIES_CONFIG), collapse = ", "))
