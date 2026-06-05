# =============================================================================
# Top correlations with half-life
# =============================================================================
# For every numeric column, compute its Spearman correlation with halflife,
# BH-correct the p-values, and plot the top-N strongest correlations.
#
# Usage (as a function — reusable):
#   source("R/load_all.R")
#   source("analysis/correlations/halflife_correlation.R")
#   df  <- build_dataset("human")
#   out <- halflife_correlation_plot(df, top_n = 20)
#   print(out$plot)
#   head(out$table)

# Or run this file top-to-bottom — the block at the bottom produces a default
# plot for human and saves it to data/outputs/plots/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(forcats)
})


#' Compute Spearman correlations of every numeric column with a response
#'
#' @param df A dataframe.
#' @param response Character. Column name of the response variable.
#' @param exclude Character vector. Regex patterns for columns to skip. Useful
#'   for dropping predictions/derived columns that would be circular.
#' @return A tibble with columns variable, correlation, p_value, p_adj_bh,
#'   abs_corr, sorted by abs_corr descending.
#' @export
correlate_with_response <- function(df,
                                    response = "halflife",
                                    exclude  = c("^saluki_prediction$",
                                                 "^prediction_difference$")) {
  if (!response %in% names(df)) stop("response column '", response, "' not found")

  numeric_df <- df |> dplyr::select(where(is.numeric))
  x_vars <- setdiff(names(numeric_df), response)

  # Apply exclusion patterns
  for (pat in exclude) x_vars <- x_vars[!grepl(pat, x_vars)]

  y <- numeric_df[[response]]

  purrr::map_df(x_vars, function(var) {
    res <- tryCatch(
      cor.test(numeric_df[[var]], y, method = "spearman", exact = FALSE),
      error   = function(e) NULL,
      warning = function(w) {
        suppressWarnings(cor.test(numeric_df[[var]], y,
                                  method = "spearman", exact = FALSE))
      }
    )
    if (is.null(res)) {
      tibble::tibble(variable = var, correlation = NA_real_, p_value = NA_real_)
    } else {
      tibble::tibble(variable = var,
                     correlation = unname(res$estimate),
                     p_value = res$p.value)
    }
  }) |>
    dplyr::mutate(
      p_adj_bh = p.adjust(p_value, method = "BH"),
      abs_corr = abs(correlation)
    ) |>
    dplyr::arrange(dplyr::desc(abs_corr))
}


#' Plot the top-N strongest correlations with half-life
#'
#' @param df A dataframe.
#' @param response Character. Column name of the response (default "halflife").
#' @param top_n Integer. How many features to show (default 20).
#' @param exclude Character vector of regex patterns to skip (passed to
#'   correlate_with_response).
#' @param formatter Function. Used to render column names on the y-axis.
#' @return A list with two elements:
#'   - `plot`: the ggplot object
#'   - `table`: the full correlation tibble (variable, correlation, p_value,
#'     p_adj_bh, abs_corr) sorted by abs_corr descending
#' @export
halflife_correlation_plot <- function(df,
                                      response = "halflife",
                                      top_n    = 20,
                                      exclude  = c("^saluki_prediction$",
                                                   "^prediction_difference$"),
                                      formatter = format_col_name) {
  corr <- correlate_with_response(df, response = response, exclude = exclude)

  top <- corr |>
    dplyr::filter(!is.na(correlation)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      q_label = paste0("q = ", format(p_adj_bh, scientific = TRUE, digits = 2)),
      c_label   = paste0(round(correlation, 2), " (", q_label, ")")
    )
  stopifnot("c_label" %in% names(top))
  p <- ggplot(top,
              aes(x = abs_corr, y = forcats::fct_reorder(variable, abs_corr))) +
    geom_col(aes(fill = correlation > 0), show.legend = FALSE) +
    geom_text(aes(label = c_label), hjust = -0.1, size = 3.5) +
    scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "firebrick")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    scale_y_discrete(labels = formatter) +
    labs(
      title    = paste0("Top ", top_n, " Spearman correlations with ",
                        formatter(response)),
      subtitle = "BH-adjusted q-values in parentheses",
      x        = "Spearman correlation coefficient",
      y        = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title       = element_text(size = 20, face = "bold", hjust = 1),
      plot.subtitle    = element_text(size = 14, hjust = 1),
      axis.title.x     = element_text(size = 18),
      axis.text        = element_text(size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, table = corr)
}


# --- Top-to-bottom run -------------------------------------------------------
# Runs only when this file is sourced/executed directly.
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  if (exists("build_dataset")) {
    df  <- build_dataset("human")
    out <- halflife_correlation_plot(df, top_n = 20)
    print(out$plot)
    
    dir.create(file.path(OUTPUT_DIR, "plots"),  showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)
    
    ggsave(file.path(OUTPUT_DIR, "plots", "halflife_correlation.jpg"),
           plot = out$plot, width = 210, height = 148, units = "mm", dpi = 300)
    write.csv(out$table,
              file.path(OUTPUT_DIR, "tables", "halflife_spearman.csv"),
              row.names = FALSE)
  }
}
