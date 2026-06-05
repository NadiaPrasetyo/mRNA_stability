# =============================================================================
# LASSO model for half-life
# =============================================================================
# Fits a LASSO regression of halflife on a user-specified bundle of features,
# with cross-validated lambda. Designed to be species-agnostic: the dataframe
# decides the species.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/models/lasso.R")
#   df   <- build_dataset("human")
#   fit  <- run_lasso(df)
#   print(fit$r_squared)
#   print(lasso_coefficient_plot(fit, top_n = 12))
# =============================================================================

suppressPackageStartupMessages({
  library(glmnet)
  library(dplyr)
  library(ggplot2)
  library(forcats)
})


# --- Default feature bundle --------------------------------------------------
# Edit this to change what goes into the model. Each entry is either a
# feature-group key (matched through FEATURE_PATTERNS) or a literal column
# name. The helper `assemble_predictors()` resolves them against a dataframe.

DEFAULT_LASSO_PREDICTORS <- list(
  groups = c("lengths",   "lengths_some",  "gc",   "nmd", "architecture",   
                  "junctions",      
              "uorfs",   "stopfree",    "skews",   "distances",      
            "codon_freqs", "aa_freqs", "nuc_ratios"),
  # groups = c("lengths",         "lengths_some",    "gc",              "nmd",             "architecture",   
  #            "mfe_scores",      "mfe_zscores",     "rnafold_scores",  "rnafold_zscores", "mfe_deltas",     
  #            "mfe_expected",    "rnalfold_scores",     "rnalfold_zscores",    "rnaup",           "junctions",      
  #            "uorfs",   "stopfree",        "skews",           "distances",      
  #            "codon_freqs", "aa_freqs", "nuc_ratios"),
  columns = c()
)


#' Resolve a predictor spec (groups + explicit columns) against a dataframe
#'
#' @param df A dataframe.
#' @param spec A list with elements `groups` (character) and `columns`
#'   (character). Missing columns are silently dropped with a message.
#' @return Character vector of columns that exist in df.
#' @export
assemble_predictors <- function(df, spec = DEFAULT_LASSO_PREDICTORS) {
  cols_from_groups <- character(0)
  for (g in spec$groups %||% character(0)) {
    cols_from_groups <- c(cols_from_groups, fg_columns(df, g))
  }
  requested <- unique(c(cols_from_groups, spec$columns %||% character(0)))
  missing   <- setdiff(spec$columns %||% character(0), names(df))
  if (length(missing) > 0) {
    message("Dropping requested columns not in df: ",
            paste(missing, collapse = ", "))
  }
  intersect(requested, names(df))
}


`%||%` <- function(a, b) if (is.null(a)) b else a


#' Fit a cross-validated LASSO model for the response
#'
#' @param df A dataframe (typically from build_dataset()).
#' @param response Column name of the response variable (default "halflife").
#' @param predictor_spec A list(groups, columns) — see DEFAULT_LASSO_PREDICTORS.
#' @param train_frac Fraction of data for training (rest is held out).
#' @param nfolds Folds for CV.
#' @param seed Integer seed for reproducibility.
#' @return A list with elements:
#'   - model       : the glmnet fit
#'   - cv          : the cv.glmnet fit
#'   - best_lambda : the chosen lambda (1se)
#'   - coefficients: tidy tibble (Feature, Coefficient)
#'   - r_squared   : held-out R-squared
#'   - rmse        : held-out RMSE
#'   - n_train, n_test
#' @export
run_lasso <- function(df,
                      response       = "halflife",
                      predictor_spec = DEFAULT_LASSO_PREDICTORS,
                      train_frac     = 0.8,
                      nfolds         = 10,
                      seed           = 12) {
  if (!response %in% names(df)) stop("response '", response, "' not found")

  predictors <- assemble_predictors(df, predictor_spec)
  if (length(predictors) == 0) stop("No predictors resolved from spec.")

  model_df <- df |>
    dplyr::select(dplyr::all_of(c(response, predictors))) |>
    na.omit()

  message("Complete cases: ", nrow(model_df),
          " / ", nrow(df), " (", length(predictors), " predictors)")

  y <- model_df[[response]]
  X <- model.matrix(~ . - 1,
                    data = model_df[, predictors, drop = FALSE])

  set.seed(seed)
  idx_train <- sample(seq_len(nrow(X)), floor(train_frac * nrow(X)))
  X_train <- X[idx_train, , drop = FALSE]
  X_test  <- X[-idx_train, , drop = FALSE]
  y_train <- y[idx_train]
  y_test  <- y[-idx_train]

  cv_fit <- glmnet::cv.glmnet(X_train, y_train, alpha = 1, nfolds = nfolds)
  best_lambda <- cv_fit$lambda.1se
  final_fit   <- glmnet::glmnet(X_train, y_train, alpha = 1,
                                lambda = best_lambda)

  coefs <- coef(final_fit, s = "lambda.1se")
  coef_df <- tibble::tibble(
    Feature     = coefs@Dimnames[[1]][coefs@i + 1],
    Coefficient = coefs@x
  )

  preds <- predict(final_fit, newx = X_test)
  r_squared <- 1 - sum((y_test - preds)^2) /
                   sum((y_test - mean(y_test))^2)
  rmse <- sqrt(mean((preds - y_test)^2))

  list(
    model        = final_fit,
    cv           = cv_fit,
    best_lambda  = best_lambda,
    coefficients = coef_df,
    r_squared    = r_squared,
    rmse         = rmse,
    n_train      = length(y_train),
    n_test       = length(y_test),
    response     = response
  )
}


#' Plot the top-N LASSO coefficients by absolute magnitude
#'
#' @param fit Result of `run_lasso()`.
#' @param top_n Number of features to display.
#' @param formatter Function to render feature names.
#' @return A ggplot object.
#' @export
lasso_coefficient_plot <- function(fit,
                                   top_n = 12,
                                   formatter = format_col_name) {
  coefs <- fit$coefficients |>
    dplyr::filter(Feature != "(Intercept)") |>
    dplyr::rowwise() |>
    dplyr::mutate(Feature = formatter(Feature),
                  label = as.character(round(Coefficient, 3))) |>
    dplyr::ungroup() |>
    dplyr::slice_max(order_by = abs(Coefficient), n = top_n)

  ggplot(coefs,
         aes(x = abs(Coefficient),
             y = forcats::fct_reorder(Feature, abs(Coefficient)))) +
    geom_col(aes(fill = Coefficient > 0), show.legend = FALSE) +
    geom_text(aes(label = label), hjust = -0.15, size = 5) +
    scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "firebrick")) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
    labs(
      title    = paste0("LASSO model for ", formatter(fit$response)),
      subtitle = "Blue = positive (stabilising), red = negative (destabilising)",
      x = "Absolute coefficient", y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title       = element_text(size = 22, face = "bold", hjust = 1),
      plot.subtitle    = element_text(size = 14, hjust = 1),
      axis.title.x     = element_text(size = 18),
      axis.text        = element_text(size = 16),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )
}


# --- Top-to-bottom run -------------------------------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  if (exists("build_dataset")) {
    df  <- build_dataset("human")
    fit <- run_lasso(df)

    cat(sprintf("Held-out R^2: %.4f\n", fit$r_squared))
    cat(sprintf("Held-out RMSE: %.4f\n", fit$rmse))

    p <- lasso_coefficient_plot(fit, top_n = 12)
    print(p)

    dir.create(file.path(OUTPUT_DIR, "plots"), showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(OUTPUT_DIR, "plots", "lasso_halflife.jpg"),
           plot = p, width = 210, height = 148, units = "mm", dpi = 300)
  }
}
