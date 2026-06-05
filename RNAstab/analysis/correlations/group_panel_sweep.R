# =============================================================================
# Feature-group correlation panel sweep
# =============================================================================
# Iterates every group in FEATURE_PATTERNS, producing a faceted scatter panel
# of group members vs the response (default halflife) plus a tibble of per-
# variable correlation statistics. Writes one plot + table per group, plus a
# combined cross-group summary table.
#
# This file bundles two functions:
#
#   feature_group_panel(df, group, ...)        — per §7.1 of PIPELINE_GUIDE,
#                                                  works on a single group.
#   feature_group_panel_sweep(df, groups, ...) — iterates groups, returns a
#                                                  named list of results.
#
# If you already have feature_group_panel deployed from §7.1, delete the copy
# below — the sweep will use whichever is in scope.
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/group_panel_sweep.R")
#   df  <- build_all()                                # all species, faceted
#   res <- feature_group_panel_sweep(df)             # all groups
#   res$rnafold_zscores$plot                          # inspect one
#
#   # single species (no species facet row):
#   df  <- build_dataset("human")
#   res <- feature_group_panel_sweep(df)
#
# Or run end-to-end:
#   Rscript analysis/correlations/group_panel_sweep.R
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


# -----------------------------------------------------------------------------
# Default-view ergonomics (analysis-layer config, not schema)
# -----------------------------------------------------------------------------
# Groups excluded from the default sweep (groups = NULL). These are real,
# distinct FEATURE_PATTERNS families — not aliases or subsets — but their high
# cardinality makes them unsuitable for an at-a-glance one-panel-per-group
# sweep (a 64-facet codon panel would be truncated by max_facets and silently
# hide members). Request them explicitly with e.g.
#   feature_group_panel_sweep(df, groups = "codon_freqs", max_facets = 64)
# This is specific to this tool's default view, so it lives here, not config.R.
DEFAULT_SWEEP_SKIP <- c("codon_freqs", "aa_freqs")


# -----------------------------------------------------------------------------
# Single-group worker (mirrors §7.1 worked example)
# -----------------------------------------------------------------------------

#' Faceted scatter panel of every column in a feature group vs response.
#'
#' @param df          Dataframe from build_dataset() or build_all().
#' @param group       Character. A key of FEATURE_PATTERNS.
#' @param response    Character. Response column (default "halflife").
#' @param method      Correlation method (default "spearman").
#' @param max_facets  Integer. If group has more members than this, only the
#'                    top-N by |correlation| are plotted; the table keeps all.
#' @param hex_bins    Integer. Number of bins per axis for geom_hex (default 50).
#' @param formatter   Display formatter (default format_col_name).
#' @return list(plot, table). Table: variable, n, correlation, p_value,
#'   p_adj_bh, abs_correlation, rank, plotted (logical).
#' @export
feature_group_panel <- function(df,
                                group,
                                response   = "halflife",
                                method     = c("spearman", "pearson",
                                               "kendall"),
                                max_facets = 25,
                                hex_bins   = 50,
                                formatter  = format_col_name) {
  
  method <- match.arg(method)
  
  # --- R5: guard ----------------------------------------------------------
  if (!response %in% names(df)) {
    stop("response '", response, "' not in df")
  }
  if (!group %in% names(FEATURE_PATTERNS)) {
    extra <- ""
    if (group %in% names(SUPERGROUPS)) {
      extra <- " (it is a supergroup; this tool is single-group only)"
    } else if (exists("GROUP_BUNDLES", inherits = TRUE) &&
               group %in% names(GROUP_BUNDLES)) {
      extra <- " (it is a bundle; this tool is single-group only)"
    }
    stop("group '", group, "' not in FEATURE_PATTERNS", extra,
         " — see R/config.R")
  }
  
  # --- R3: use fg_columns to enumerate ------------------------------------
  cols <- fg_columns(df, group)
  if (length(cols) == 0) {
    stop("group '", group, "' resolved to 0 columns in this df")
  }
  
  # --- Long-form ----------------------------------------------------------
  long <- df |>
    dplyr::select(dplyr::any_of(c("species", response)),
                  dplyr::all_of(cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(cols),
                        names_to  = "variable",
                        values_to = "value") |>
    dplyr::filter(!is.na(value), !is.na(.data[[response]]))
  
  if (nrow(long) == 0) {
    stop("group '", group, "' has no non-NA pairs with response")
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
  
  # --- Select what to plot ------------------------------------------------
  truncated <- nrow(summary_tbl) > max_facets
  plot_vars <- if (truncated) {
    summary_tbl$variable[seq_len(max_facets)]
  } else {
    summary_tbl$variable
  }
  
  summary_tbl <- summary_tbl |>
    dplyr::mutate(plotted = variable %in% plot_vars)
  
  long_plot <- long |> dplyr::filter(variable %in% plot_vars)
  
  # --- R4: facet labels via formatter -------------------------------------
  # Append rho + q to the strip label so each panel is self-documenting.
  annot <- summary_tbl |>
    dplyr::filter(plotted) |>
    dplyr::mutate(strip = sprintf(
      "%s\nrho = %+.2f, q = %.2g",
      formatter(variable), correlation, p_adj_bh
    ))
  strip_labels <- setNames(annot$strip, annot$variable)
  
  # Order facets by signed correlation (most stabilising first)
  long_plot$variable <- factor(
    long_plot$variable,
    levels = summary_tbl |>
      dplyr::filter(plotted) |>
      dplyr::arrange(dplyr::desc(correlation)) |>
      dplyr::pull(variable)
  )
  
  # --- Plot ---------------------------------------------------------------
  # Hex density for the point cloud (count on a log10 viridis scale), with the
  # lm trend drawn on top. When multiple species are present we facet the
  # density surface by species (one row per species) rather than pooling, so
  # each species gets its own count colourbar-comparable surface.
  p <- ggplot2::ggplot(long_plot,
                       ggplot2::aes(x = value, y = .data[[response]])) +
    ggplot2::geom_hex(bins = hex_bins) +
    ggplot2::scale_fill_viridis_c(trans = "log10", name = "count") +
    ggplot2::geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE,
      colour = "firebrick", linewidth = 0.6
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
      title    = sprintf("%s vs %s — group: %s",
                         formatter(response),
                         tools::toTitleCase(method),
                         group),
      subtitle = if (truncated) {
        sprintf("%d variables, showing top %d by |%s|; full ranking in table",
                nrow(summary_tbl), max_facets, method)
      } else {
        sprintf("%d variables, all shown; %s with BH-adjusted q",
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
# Sweep wrapper — iterate FEATURE_PATTERNS
# -----------------------------------------------------------------------------

#' Run feature_group_panel over every group in FEATURE_PATTERNS.
#'
#' One panel per group. Groups must be FEATURE_PATTERNS keys — this tool shows
#' each schema family on its own, so it does not accept supergroups or bundles
#' (those are intent-layer objects for multi-group / refined selections; use
#' feature_correlation_dotplot() or feature_response_scatter() for those).
#'
#' Default (groups = NULL): every family minus DEFAULT_SWEEP_SKIP (high-
#' cardinality groups like codon_freqs/aa_freqs). Groups that resolve to zero
#' columns are skipped with a message rather than erroring.
#'
#' @param df         Dataframe from build_dataset() or build_all().
#' @param groups     Character vector of FEATURE_PATTERNS keys. Default: all
#'                   real families minus DEFAULT_SWEEP_SKIP.
#' @param response   Character. Response column (default "halflife").
#' @param method     Correlation method (default "spearman").
#' @param max_facets Integer cap on facets per panel (default 25).
#' @param hex_bins   Integer bins per axis for geom_hex (default 50).
#' @param formatter  Display formatter (default format_col_name).
#' @return Named list keyed by group; each element is list(plot, table).
#' @export
feature_group_panel_sweep <- function(df,
                                      groups     = NULL,
                                      response   = "halflife",
                                      method     = "spearman",
                                      max_facets = 25,
                                      hex_bins   = 50,
                                      formatter  = format_col_name) {
  
  if (!response %in% names(df)) stop("response '", response, "' not in df")
  
  if (is.null(groups)) {
    # Every real FEATURE_PATTERNS family, minus high-cardinality groups that
    # don't suit an at-a-glance sweep (see DEFAULT_SWEEP_SKIP). Aliases/subsets
    # no longer exist in the schema, so cardinality is the only skip reason.
    groups <- setdiff(expand_groups(NULL), DEFAULT_SWEEP_SKIP)
  }
  
  results <- list()
  
  for (g in groups) {
    cols <- fg_columns(df, g)
    if (length(cols) == 0) {
      message("Skipping group '", g, "': 0 matching columns in df")
      next
    }
    message("Building panel for '", g, "' (", length(cols), " columns)...")
    
    # Defensive: if a single group fails (e.g. all NA after filter), record
    # the skip and carry on — never let one group abort the sweep.
    res <- tryCatch(
      feature_group_panel(df, g,
                          response   = response,
                          method     = method,
                          max_facets = max_facets,
                          hex_bins   = hex_bins,
                          formatter  = formatter),
      error = function(e) {
        message("  failed: ", conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(res)) results[[g]] <- res
  }
  
  results
}


# -----------------------------------------------------------------------------
# Helper: pick A4-ish dimensions based on facet count
# -----------------------------------------------------------------------------

# Picks output (width, height) in mm based on how many facets the plot will
# contain. Keeps small groups compact, large groups readable.
#
# n_species > 1 signals a facet_grid(species ~ variable) layout: every variable
# is its own column and every species its own row, so width tracks n_facets
# directly (capped) and height tracks n_species. n_species == 1 keeps the
# original facet_wrap sizing.
panel_dimensions <- function(n_facets, n_species = 1, ncol_grid = NULL) {
  if (n_species > 1) {
    ncol_grid <- min(n_facets, 8L)   # cap columns; wide panels stay legible
    n_rows    <- n_species
    width  <- 80 + ncol_grid * 50
    height <- 80 + n_rows    * 55
    return(list(width = width, height = height))
  }
  
  if (is.null(ncol_grid)) {
    ncol_grid <- if      (n_facets <= 3)  n_facets
    else if (n_facets <= 6)  3
    else if (n_facets <= 12) 4
    else                     5
  }
  n_rows <- ceiling(n_facets / ncol_grid)
  width  <- 80  + ncol_grid * 50
  height <- 80  + n_rows    * 50
  list(width = width, height = height)
}


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  
  df  <- build_all()
  res <- feature_group_panel_sweep(df)
  
  if (length(res) == 0) stop("Sweep produced no results — check the dataset")
  
  n_species <- if ("species" %in% names(df)) {
    dplyr::n_distinct(df$species)
  } else 1L
  
  # R8: outputs under OUTPUT_DIR
  dir.create(file.path(OUTPUT_DIR, "plots"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"),
             showWarnings = FALSE, recursive = TRUE)
  
  for (g in names(res)) {
    n_plotted <- sum(res[[g]]$table$plotted)
    dims      <- panel_dimensions(n_plotted, n_species = n_species)
    
    ggplot2::ggsave(
      file.path(OUTPUT_DIR, "plots",
                paste0("group_panel_", g, ".jpg")),
      plot   = res[[g]]$plot,
      width  = dims$width, height = dims$height,
      units = "mm", dpi = 300, limitsize = FALSE
    )
    write.csv(
      res[[g]]$table,
      file.path(OUTPUT_DIR, "tables",
                paste0("group_panel_", g, ".csv")),
      row.names = FALSE
    )
  }
  
  # Cross-group summary: one row per (group, variable), with the group label
  # attached. Useful for "which group carries the most signal" questions
  # without re-running anything.
  combined <- purrr::imap_dfr(res, function(r, g) {
    r$table |> dplyr::mutate(group = g, .before = 1)
  })
  write.csv(
    combined,
    file.path(OUTPUT_DIR, "tables", "group_panel_combined_summary.csv"),
    row.names = FALSE
  )
  
  message("\nSweep complete:")
  message("  ", length(res), " group panels in ",
          file.path(OUTPUT_DIR, "plots"))
  message("  ", length(res), " per-group tables + 1 combined summary in ",
          file.path(OUTPUT_DIR, "tables"))
  message("\nQuick eyeball — top variable per group by |rho|:")
  combined |>
    dplyr::filter(rank == 1) |>
    dplyr::select(group, variable, correlation, n, p_adj_bh) |>
    dplyr::arrange(dplyr::desc(abs(correlation))) |>
    print(n = Inf)
}