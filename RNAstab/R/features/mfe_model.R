# =============================================================================
# MFE thermodynamic model
# =============================================================================
# Parameters from in-house model
#
# Expected MFE = (a * (GC/100)^b + c) * length + d,  capped at 0.
# =============================================================================


#' Calculate expected MFE for given GC content and sequence length
#'
#' Vectorised. GC is expected to be in percent (0-100).
#'
#' @param gc Numeric vector of GC content in percent.
#' @param length Numeric vector of sequence lengths (nt).
#' @return Numeric vector of expected MFE values. Capped at 0 (folding free
#'   energy cannot be positive).
#' @export
calculate_mfe_expected <- function(gc, length) {
  a <- -0.8403211
  b <-  2.3521348
  c <- -0.1534111
  d <- 13.5600933
  gc_frac <- gc / 100
  expected <- (a * gc_frac^b + c) * length + d
  pmin(expected, 0)
}


#' Calculate the thermodynamic delta between observed and expected MFE
#'
#' @param mfe Observed MFE score.
#' @param mfe_expected Expected MFE for the same (GC, length).
#' @return mfe - mfe_expected. Proportional to the log-odds of the folded
#'   state given sequence composition.
#' @export
calculate_mfe_delta <- function(mfe, mfe_expected) mfe - mfe_expected
