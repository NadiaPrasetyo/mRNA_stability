# =============================================================================
# Cache I/O
# =============================================================================
# Save and load built dataset snapshots as .rds files. Filenames include a
# version integer (CACHE_VERSION in config.R) so stale caches are naturally
# replaced when feature-engineering logic changes.
# =============================================================================


#' Save a built dataset to the cache.
#'
#' @param df A dataframe.
#' @param species Character, one of names(SPECIES_CONFIG).
#' @return The path the file was written to, invisibly.
#' @export
save_snapshot <- function(df, species) {
  dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  path <- cache_path(species)
  saveRDS(df, path)
  message(sprintf("Cache saved: %s (%d rows, %d columns)",
                  path, nrow(df), ncol(df)))
  invisible(path)
}


#' Load a built dataset from the cache.
#'
#' @param species Character, one of names(SPECIES_CONFIG).
#' @return The cached dataframe, or NULL if no cache exists.
#' @export
load_snapshot <- function(species) {
  path <- cache_path(species)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}


#' Remove the cached dataset for a species.
#'
#' @param species Character, one of names(SPECIES_CONFIG).
#' @return TRUE if a cache was deleted, FALSE if no cache existed.
#' @export
clear_snapshot <- function(species) {
  path <- cache_path(species)
  if (file.exists(path)) {
    file.remove(path)
    message("Cache cleared: ", path)
    return(TRUE)
  }
  FALSE
}
