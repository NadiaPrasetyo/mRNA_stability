# =============================================================================
# Feature × response hex-density panels (selection-driven)
# =============================================================================
# A middle ground between scatter_plot.R (one x vs one y, per-point colour) and
# group_panel_sweep.R (one group, one facet per column). This script:
#
#   * takes a SELECTION (groups / supergroups / bundles / individual columns),
#     not a single FEATURE_PATTERNS key, via the selection layer;
#   * draws one hex-density panel per selected feature — feature on x, the
#     response (default `halflife`) on y — with NO third colour dimension
#     (fill = bin count, log10 viridis, exactly like group_panel);
#   * overlays BOTH an lm (firebrick) and a loess (steel-blue) trend per panel;
#   * annotates each strip with the per-feature rho and BH-adjusted q.
#
# Unlike feature_group_panel(), this accepts the full selection vocabulary, so
# you can hex-panel an ad-hoc set:
#   feature_response_hex_panels(df, columns = c("length_cds", "gc_content_cds"))
#   feature_response_hex_panels(df, groups  = "nmd_reported")          # a bundle
#   feature_response_hex_panels(df, groups  = "structure",
#                               drop = list(probing = "gini_nucleoplasm_cds"))
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_response_hex_panels.R")
#   df  <- build_dataset("human")
#   out <- feature_response_hex_panels(df, groups = "intrinsic")
#   print(out$plot); head(out$table)
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(forcats)
  library(hexbin)   # required by ggplot2::geom_hex()
})


#' Hex-density panels of each selected feature vs a response.
#'
#' Selection is resolved through the project selection layer, so `groups`
#' accepts FEATURE_PATTERNS keys, SUPERGROUPS names, and GROUP_BUNDLES names;
#' `columns` accepts literal column names. The plotted feature set is the union
#' of both, after pick/drop refinement, with the response and any derived
#' predictions excluded (R10).
#'
#' @param df         Dataframe from build_dataset() or build_all().
#' @param response   Character. Response column on the y-axis (default
#'                   "halflife"). Any numeric column is permitted.
#' @param groups     Character vector of group / supergroup / bundle names, or
#'                   NULL. NULL with `columns = NULL` selects every group.
#' @param columns    Character vector of literal column names to include
#'                   alongside whatever `groups` resolves to. Union semantics.
#' @param pick       Named list: group key -> columns to keep (allow-list).
#' @param drop       Named list: group key -> columns to remove.
#' @param method     Correlation method for the strip annotation (default
#'                   "spearman").
#' @param max_facets Integer. If the selection resolves to more features than
#'                   this, only the top-N by |correlation| are plotted; the
#'                   table keeps all. Default 25.
#' @param hex_bins   Integer bins per axis for geom_hex (default 50).
#' @param exclude    Regex patterns to drop from the feature set (R10). The
#'                   response is always excluded regardless.
#' @param formatter  Display formatter (default format_col_name).
#' @return list(plot, table). Table: variable, n, correlation, p_value,
#'   p_adj_bh, abs_correlation, rank, plotted (logical).
#' @export
feature_response_hex_panels <- function(df,
                                        response   = "halflife",
                                        groups     = NULL,
                                        columns    = NULL,
                                        pick       = list(),
                                        drop       = list(),
                                        method     = c("spearman", "pearson",
                                                       "kendall"),
                                        max_facets = 25,
                                        hex_bins   = 50,
                                        exclude    = c("^saluki_prediction$",
                                                       "^prediction_difference$"),
                                        formatter  = format_col_name) {
  
  method <- match.arg(method)
  
  # --- R5: guards ---------------------------------------------------------
  if (!response %in% names(df)) {
    stop("response '", response, "' not in df")
  }
  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }
  
  # --- Resolve the feature set (selection layer ∪ literal columns) --------
  # When the caller names nothing at all, default groups to every group; but if
  # they passed only `columns`, respect that and don't drag in all groups.
  sel_cols <- if (is.null(groups) && !is.null(columns)) {
    character()
  } else {
    select_features(df, groups = groups, pick = pick, drop = drop)
  }
  
  lit_cols <- character()
  if (!is.null(columns)) {
    missing_lit <- setdiff(columns, names(df))
    if (length(missing_lit) > 0) {
      message("feature_response_hex_panels: columns not in df (ignored): ",
              paste(missing_lit, collapse = ", "))
    }
    lit_cols <- intersect(columns, names(df))
  }
  
  candidates <- union(sel_cols, lit_cols)
  
  # R10 + response: drop derived predictions and the response itself.
  for (rgx in exclude) {
    candidates <- candidates[!grepl(rgx, candidates)]
  }
  candidates <- setdiff(candidates, response)
  
  # Keep only numeric features (hex needs a numeric x).
  candidates <- candidates[vapply(candidates,
                                  function(c) is.numeric(df[[c]]),
                                  logical(1))]
  
  if (length(candidates) == 0) {
    stop("No numeric feature columns after selection/exclusion — check ",
         "`groups`, `columns`, and `exclude`")
  }
  
  # --- Long-form ----------------------------------------------------------
  long <- df |>
    dplyr::select(dplyr::any_of(c("species", response)),
                  dplyr::all_of(candidates)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(candidates),
                        names_to  = "variable",
                        values_to = "value") |>
    dplyr::filter(!is.na(value), !is.na(.data[[response]]))
  
  if (nrow(long) == 0) {
    stop("No non-NA pairs between the selected features and '", response, "'")
  }
  
  has_species <- "species" %in% names(long) &&
    length(unique(long$species)) > 1
  
  # --- R9: per-variable summary table -------------------------------------
  summary_tbl <- long |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      n           = dplyr::n(),
      correlation = suppressWarnings(stats::cor(
        value, .data[[response]],
        method = method, use = "pairwise.complete.obs"
      )),
      p_value     = suppressWarnings(stats::cor.test(
        value, .data[[response]],
        method = method, exact = FALSE
      )$p.value),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      p_adj_bh        = stats::p.adjust(p_value, method = "BH"),
      abs_correlation = abs(correlation),
      rank            = dplyr::min_rank(dplyr::desc(abs_correlation))
    ) |>
    dplyr::arrange(rank)
  
  # --- Select what to plot (table keeps all) ------------------------------
  truncated <- nrow(summary_tbl) > max_facets
  plot_vars <- if (truncated) {
    summary_tbl$variable[seq_len(max_facets)]
  } else {
    summary_tbl$variable
  }
  summary_tbl <- summary_tbl |>
    dplyr::mutate(plotted = variable %in% plot_vars)
  
  long_plot <- long |> dplyr::filter(variable %in% plot_vars)
  
  # --- R4: strip labels via formatter, annotated with rho + q -------------
  annot <- summary_tbl |>
    dplyr::filter(plotted) |>
    dplyr::mutate(strip = sprintf(
      "%s\nrho = %+.2f, q = %.2g",
      formatter(variable), correlation, p_adj_bh
    ))
  strip_labels <- stats::setNames(annot$strip, annot$variable)
  
  # Order facets by signed correlation (most stabilising first).
  long_plot$variable <- factor(
    long_plot$variable,
    levels = summary_tbl |>
      dplyr::filter(plotted) |>
      dplyr::arrange(dplyr::desc(correlation)) |>
      dplyr::pull(variable)
  )
  
  # --- Plot ---------------------------------------------------------------
  # Hex density (count on a log10 viridis scale) with BOTH trend lines: lm
  # (firebrick) and loess (steel-blue). No third colour dimension — fill is
  # purely the bin count, matching group_panel.
  p <- ggplot2::ggplot(long_plot,
                       ggplot2::aes(x = value, y = .data[[response]])) +
    ggplot2::geom_hex(bins = hex_bins) +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "count") +
    ggplot2::geom_smooth(
      ggplot2::aes(colour = "lm"),
      method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6
    ) +
    ggplot2::geom_smooth(
      ggplot2::aes(colour = "loess"),
      method = "loess", formula = y ~ x, se = FALSE,
      linewidth = 0.6, span = 0.75
    ) +
    ggplot2::scale_colour_manual(
      name   = "trend",
      values = c(lm = "firebrick", loess = "steelblue"),
      labels = c(lm = "Linear (lm)", loess = "LOESS")
    )
  
  if (has_species) {
    p <- p + ggplot2::facet_grid(
      species ~ variable, scales = "free_x",
      labeller = ggplot2::labeller(
        variable = strip_labels,
        species  = ggplot2::as_labeller(formatter)
      )
    )
  } else {
    p <- p + ggplot2::facet_wrap(
      ~ variable, scales = "free_x",
      labeller = ggplot2::labeller(variable = strip_labels)
    )
  }
  
  p <- p +
    ggplot2::labs(
      title    = sprintf("%s vs selected features — %s",
                         formatter(response),
                         tools::toTitleCase(method)),
      subtitle = if (truncated) {
        sprintf("%d features selected, showing top %d by |%s|; full ranking in table",
                nrow(summary_tbl), max_facets, method)
      } else {
        sprintf("%d features, all shown; %s with BH-adjusted q",
                nrow(summary_tbl), tools::toTitleCase(method))
      },
      x = NULL,
      y = formatter(response)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 10),
      strip.text       = ggplot2::element_text(size = 8),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right"
    )
  
  list(plot = p, table = summary_tbl)
}


# -----------------------------------------------------------------------------
# Helper: A4-ish dimensions based on facet count (mirrors group_panel_sweep)
# -----------------------------------------------------------------------------
panel_dimensions <- function(n_facets, n_species = 1, ncol_grid = NULL) {
  if (n_species > 1) {
    ncol_grid <- min(n_facets, 8L)
    width  <- 80 + ncol_grid * 50
    height <- 80 + n_species  * 55
    return(list(width = width, height = height))
  }
  if (is.null(ncol_grid)) {
    ncol_grid <- if      (n_facets <= 3)  n_facets
    else if (n_facets <= 6)  3
    else if (n_facets <= 12) 4
    else                     5
  }
  n_rows <- ceiling(n_facets / ncol_grid)
  width  <- 80 + ncol_grid * 50
  height <- 80 + n_rows    * 50
  list(width = width, height = height)
}


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------
if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  
  df <- build_dataset("human")
  
  dir.create(file.path(OUTPUT_DIR, "plots"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"),
             showWarnings = FALSE, recursive = TRUE)
  
  n_species <- if ("species" %in% names(df)) dplyr::n_distinct(df$species) else 1L
  
  jobs <- list(
    list(suffix = "intrinsic", groups = "intrinsic", columns = NULL),
    list(suffix = "structure", groups = "structure", columns = NULL)
  )
  
  for (job in jobs) {
    message("\nHex panels for selection: ", job$suffix)
    out <- feature_response_hex_panels(
      df,
      groups  = job$groups,
      columns = job$columns
    )
    print(out$plot)
    
    n_plotted <- sum(out$table$plotted)
    dims      <- panel_dimensions(n_plotted, n_species = n_species)
    
    ggplot2::ggsave(
      file.path(OUTPUT_DIR, "plots",
                paste0("feature_response_hex_panels_", job$suffix, ".jpg")),
      plot = out$plot,
      width = dims$width, height = dims$height,
      units = "mm", dpi = 300, limitsize = FALSE
    )
    write.csv(
      out$table,
      file.path(OUTPUT_DIR, "tables",
                paste0("feature_response_hex_panels_", job$suffix, ".csv")),
      row.names = FALSE
    )
  }
  
  message("\nFeature response hex panels complete:")
  message("  ", file.path(OUTPUT_DIR, "plots"),
          "/feature_response_hex_panels_*.jpg")
}