# =============================================================================
# Normalisation functions
# =============================================================================
# Pure numeric utilities for standardising vectors.
# =============================================================================


#' Z-score normalisation
#'
#' Standardises a numeric vector so that it has mean 0 and SD 1.
#'
#' @param x A numeric vector.
#' @return Numeric vector of the same length as x. If x has fewer than 2
#'   non-NA values the result is all NA. If all non-NA values are equal
#'   the result is all 0.
#' @export
z_score_normalize <- function(x) {
  if (!is.numeric(x)) stop("Input must be a numeric vector.")

  mean_val <- mean(x, na.rm = TRUE)
  sd_val   <- sd(x, na.rm = TRUE)

  if (is.na(sd_val)) return(rep(NA_real_, length(x)))
  if (sd_val == 0)   return(rep(0, length(x)))

  (x - mean_val) / sd_val
}


#' Min-max normalisation
#'
#' Scales a numeric vector to [0, 1].
#'
#' @param x A numeric vector.
#' @return Numeric vector of the same length as x. If all non-NA values are
#'   equal the result is all 0 (NAs preserved).
#' @export
min_max_normalize <- function(x) {
  if (!is.numeric(x)) stop("Input must be a numeric vector.")

  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)

  if (min_val == max_val) {
    result <- rep(0, length(x))
    result[is.na(x)] <- NA_real_
    return(result)
  }

  (x - min_val) / (max_val - min_val)
}


#' Apply a normalisation function to every numeric column of a dataframe
#'
#' Non-numeric columns are preserved unchanged. Columns that end up entirely
#' NA after normalisation are dropped.
#'
#' @param df A dataframe.
#' @param method Either "zscore" or "minmax".
#' @return A dataframe with the same rows. Numeric columns are normalised.
#' @export
normalize_numeric <- function(df, method = c("zscore", "minmax")) {
  method <- match.arg(method)
  fn <- switch(method, zscore = z_score_normalize, minmax = min_max_normalize)

  numeric_cols <- names(df)[sapply(df, is.numeric)]
  df[numeric_cols] <- lapply(df[numeric_cols], fn)

  # Drop columns that became entirely NA
  keep <- !sapply(df, function(x) all(is.na(x)))
  df[, keep, drop = FALSE]
}
