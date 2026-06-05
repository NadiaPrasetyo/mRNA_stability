# =============================================================================
# Top-N response correlations
# =============================================================================
# Horizontal bar chart of the top-N strongest correlates of a response
# variable (default `halflife`), ranked by |correlation|. Default Spearman.
# Auto-facets by species if more than one is present.
#
# Excludes derived predictions per R10 (PIPELINE_GUIDE §4). Default exclusion:
# saluki_prediction, prediction_difference. Override via `exclude` argument
# when running against a derived response (e.g. prediction_difference itself).
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/top_correlations.R")
#
#   # Default: top-20 predictors of halflife
#   df  <- build_dataset("human")
#   out <- top_n_response_correlations(df)
#   print(out$plot)
#
#   # Bigger panel, Pearson, custom exclusion set:
#   out <- top_n_response_correlations(
#     df, top_n = 40, method = "pearson",
#     exclude = c("^saluki_prediction$", "^prediction_difference$", "^cai$")
#   )
#
#   # Repurpose for Saluki residual analysis:
#   out <- top_n_response_correlations(
#     df,
#     response = "prediction_difference",
#     exclude  = c("^saluki_prediction$", "^halflife$")
#   )
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(purrr)
  library(tibble)
})


#' Top-N predictors of a response variable, ranked by absolute correlation.
#'
#' Iterates over every numeric column in df except identifiers, the response
#' itself, and any column matching the `exclude` regex set. Computes the
#' requested correlation per species and returns the top-N by |r|.
#'
#' Facets by species when the input contains more than one (Rule R6).
#'
#' @param df         Dataframe from build_dataset() or build_all().
#' @param top_n      Integer. Number of top features to display.
#' @param response   Character. Response column (default "halflife").
#' @param method     Character. "spearman" (default), "pearson", or "kendall".
#' @param exclude    Character vector of regex patterns for columns to skip.
#'                   Default is the R10 derived-prediction exclusion set.
#' @param min_n      Integer. Minimum non-NA pairs required to compute r.
#' @param formatter  Function. Display formatter (default format_col_name).
#' @return list(plot, table). Table has columns: species, variable, n,
#'   correlation, p_value, p_adj_bh, abs_correlation, rank.
#' @export
top_n_response_correlations <- function(df,
                                        top_n     = 20,
                                        response  = "halflife",
                                        method    = c("spearman", "pearson",
                                                      "kendall"),
                                        exclude   = c("^saluki_prediction$",
                                                      "^prediction_difference$"),
                                        min_n     = 20,
                                        formatter = format_col_name) {

  method <- match.arg(method)

  # --- R5: guard ------------------------------------------------------------
  if (!response  %in% names(df)) stop("response '", response, "' not in df")
  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }

  id_cols     <- c("transcript_id", "gene_id", "gene_name", "species")
  has_species <- length(unique(df$species)) > 1

  # --- Per-species correlation worker --------------------------------------
  compute_one <- function(sub) {

    # Numeric columns only, drop identifiers and the response
    numeric_cols <- names(sub)[vapply(sub, is.numeric, logical(1))]
    candidate    <- setdiff(numeric_cols, c(id_cols, response))

    # R10: drop derived predictions
    if (length(exclude) > 0) {
      excluded_mask <- vapply(candidate, function(co) {
        any(vapply(exclude, function(rgx) grepl(rgx, co), logical(1)))
      }, logical(1))
      candidate <- candidate[!excluded_mask]
    }

    # Drop columns that are constant or near-empty in this species
    candidate <- candidate[vapply(candidate, function(co) {
      v <- sub[[co]]
      sum(!is.na(v)) >= min_n &&
        length(unique(stats::na.omit(v))) > 1
    }, logical(1))]

    if (length(candidate) == 0) {
      return(tibble::tibble(
        variable    = character(),
        n           = integer(),
        correlation = numeric(),
        p_value     = numeric()
      ))
    }

    r_vec <- sub[[response]]
    purrr::map_dfr(candidate, function(co) {
      v  <- sub[[co]]
      ok <- !is.na(v) & !is.na(r_vec)
      if (sum(ok) < min_n) return(tibble::tibble())

      ct <- suppressWarnings(
        stats::cor.test(v[ok], r_vec[ok], method = method, exact = FALSE)
      )
      tibble::tibble(
        variable    = co,
        n           = sum(ok),
        correlation = unname(ct$estimate),
        p_value     = ct$p.value
      )
    })
  }

  # --- Run per species -----------------------------------------------------
  if (has_species) {
    result <- df |>
      dplyr::group_by(species) |>
      dplyr::group_modify(~ compute_one(.x)) |>
      dplyr::ungroup()
  } else {
    result <- compute_one(df)
    if (nrow(result) > 0) result$species <- unique(df$species)
    result <- dplyr::select(result, dplyr::any_of("species"),
                            dplyr::everything())
  }

  if (nrow(result) == 0) {
    stop("No correlations could be computed — check response and exclusions")
  }

  # --- BH adjustment + ranking per species ---------------------------------
  result <- result |>
    dplyr::group_by(species) |>
    dplyr::mutate(
      p_adj_bh        = stats::p.adjust(p_value, method = "BH"),
      abs_correlation = abs(correlation),
      rank            = dplyr::min_rank(dplyr::desc(abs_correlation))
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(species, rank)

  # --- Slice top-N for the plot --------------------------------------------
  top <- result |>
    dplyr::filter(rank <= top_n) |>
    dplyr::mutate(
      display = formatter(variable),
      sign    = ifelse(correlation > 0, "Stabilising (+)", "Destabilising (-)")
    )

  # --- R4: y-axis order by signed correlation ------------------------------
  # Make display a factor ordered by signed correlation within species.
  # When faceting with free_y, ggplot will order each facet separately if we
  # use tidytext::reorder_within, but for now order globally.
  top <- top |>
    dplyr::mutate(
      display = factor(display),
      display = forcats::fct_drop(display),
      display = forcats::fct_reorder(display, correlation, .fun = mean,
                                     .na_rm = TRUE)
    )

  method_label <- c(spearman = "Spearman",
                    pearson  = "Pearson",
                    kendall  = "Kendall")[method]

  excl_pretty <- paste(gsub("[\\^\\$]", "", exclude), collapse = ", ")

  p <- ggplot2::ggplot(
    top,
    ggplot2::aes(x = correlation, y = display, fill = sign)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%+.2f", correlation),
        hjust = ifelse(correlation > 0, -0.15, 1.15)
      ),
      size = 3
    ) +
    ggplot2::scale_fill_manual(
      values = c("Stabilising (+)"   = "#2166ac",
                 "Destabilising (-)" = "#b2182b"),
      name   = "Direction"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.18, 0.18))
    ) +
    ggplot2::labs(
      title    = sprintf("Top %d predictors of %s by |%s|",
                         top_n, formatter(response), method_label),
      subtitle = if (nzchar(excl_pretty)) {
        sprintf("%s correlation; derived predictions excluded (R10): %s",
                method_label, excl_pretty)
      } else {
        sprintf("%s correlation", method_label)
      },
      x = sprintf("%s correlation with %s",
                  method_label, formatter(response)),
      y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 10),
      axis.text.y        = ggplot2::element_text(size = 10),
      legend.position    = "bottom",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (has_species) {
    p <- p + ggplot2::facet_wrap(~ species, scales = "free_y")
  }

  # --- R9: return full ranked table, not just top-N ------------------------
  # Full table is more useful for downstream consumers; the plot's subset is
  # identifiable via `rank <= top_n`.
  list(plot = p, table = result)
}


# -----------------------------------------------------------------------------
# Top-to-bottom runner
# -----------------------------------------------------------------------------

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  df <- build_dataset("human")

  out <- top_n_response_correlations(df, top_n = 20)
  print(out$plot)

  # R8: outputs under OUTPUT_DIR
  dir.create(file.path(OUTPUT_DIR, "plots"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"),
             showWarnings = FALSE, recursive = TRUE)

  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "plots", "top_correlations_halflife.jpg"),
    plot   = out$plot,
    width  = 210, height = 240, units = "mm", dpi = 300
  )
  write.csv(
    out$table,
    file.path(OUTPUT_DIR, "tables", "top_correlations_halflife.csv"),
    row.names = FALSE
  )

  message("\nTop predictors of halflife written:")
  message("  ", file.path(OUTPUT_DIR, "plots",  "top_correlations_halflife.jpg"))
  message("  ", file.path(OUTPUT_DIR, "tables", "top_correlations_halflife.csv"))
}
