# analysis/cross_species/cross_species_probing_concordance.R
# -----------------------------------------------------------------------------
# Cross-species probing concordance, per-species baseline, and ortholog
# Δ-Δ analyses for transcript-level probing metrics.
#
# Probing source: icSHAPE-derived Gini scores (gini_<compartment>_<region>
# over compartments {nucleoplasm, cytoplasm} and regions {mrna, 5utr, cds,
# 3utr}), selected via the `probing` feature group. Metric selection is
# group-driven (fg_columns(df, "probing")), so this analysis tracks whatever
# the probing group currently resolves to without code changes.
#
# Paul's question: do probing metrics that predict halflife in humans also
# do so in mouse? Decomposed into three plots, each returning list(plot, table):
#
#   1. ortholog_concordance_plot       — is the metric value itself conserved
#                                        between orthologs? (human X vs mouse X)
#   2. probing_halflife_species_plot   — does each metric predict halflife
#                                        independently in each species?
#   3. ortholog_delta_delta_plot       — does Δprobing explain Δhalflife within
#                                        ortholog pairs? (Δ = human − mouse)
#
# Ortholog pairing: by gene_name, unique within each species in all_species
# (one transcript per gene), so the human/mouse inner join is 1:1. Note the
# gini probing columns join to the main dataset on gene_id upstream (older
# annotation); cross-species orthology here is still keyed on gene_name.
# -----------------------------------------------------------------------------

# source("R/load_all.R")


# =============================================================================
# Internal helpers
# =============================================================================

#' Pair human and mouse rows by gene_name and pivot to a per-metric tidy frame.
#'
#' Δ convention: human − mouse, both for the metric (dx) and halflife (dhl).
#'
#' @param df       all_species dataframe (output of build_all()).
#' @param metrics  Character vector of metric column names present in df.
#' @return tibble with columns: gene_name, metric, x_human, x_mouse,
#'         halflife_human, halflife_mouse, dx, dhl.
pair_orthologs <- function(df, metrics) {
  # R5: guard required columns.
  required <- c("species", "gene_name", "halflife")
  missing_required <- setdiff(required, names(df))
  if (length(missing_required)) {
    stop("required columns missing from df: ",
         paste(missing_required, collapse = ", "))
  }
  missing_metrics <- setdiff(metrics, names(df))
  if (length(missing_metrics)) {
    stop("metrics not found in df: ",
         paste(missing_metrics, collapse = ", "))
  }
  if (!all(c("human", "mouse") %in% unique(df$species))) {
    stop("df must contain both 'human' and 'mouse' in the species column")
  }
  
  base_cols <- c("gene_name", "halflife", metrics)
  hum <- df |>
    dplyr::filter(species == "human", !is.na(gene_name)) |>
    dplyr::select(dplyr::all_of(base_cols))
  mou <- df |>
    dplyr::filter(species == "mouse", !is.na(gene_name)) |>
    dplyr::select(dplyr::all_of(base_cols))
  
  paired <- dplyr::inner_join(hum, mou, by = "gene_name",
                              suffix = c("_human", "_mouse"))
  
  # Long form: one row per (gene_name, metric).
  purrr::map_dfr(metrics, function(m) {
    xh <- paired[[paste0(m, "_human")]]
    xm <- paired[[paste0(m, "_mouse")]]
    tibble::tibble(
      gene_name      = paired$gene_name,
      metric         = m,
      x_human        = xh,
      x_mouse        = xm,
      halflife_human = paired$halflife_human,
      halflife_mouse = paired$halflife_mouse,
      dx             = xh - xm,
      dhl            = paired$halflife_human - paired$halflife_mouse
    )
  })
}


#' Spearman summary per group, NA-safe.
#'
#' @param df    Long data frame.
#' @param x,y   Column names (character) to correlate.
#' @param by    Character vector of grouping columns.
#' @param min_n Minimum complete pairs required to compute rho.
#' @return tibble with one row per group: <by cols>, n, rho, p_value.
spearman_summary <- function(df, x, y, by, min_n = 10) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(
      .vals = list({
        xv <- .data[[x]]
        yv <- .data[[y]]
        ok <- !is.na(xv) & !is.na(yv)
        n  <- sum(ok)
        if (n < min_n) {
          list(n = n, rho = NA_real_, p_value = NA_real_)
        } else {
          ct <- suppressWarnings(stats::cor.test(
            xv[ok], yv[ok], method = "spearman", exact = FALSE
          ))
          list(n = n, rho = unname(ct$estimate), p_value = ct$p.value)
        }
      }),
      .groups = "drop"
    ) |>
    tidyr::unnest_wider(.vals)
}


#' Format a per-panel annotation label.
annotation_label <- function(rho, p_value, n, slope = NULL) {
  is_na_rho <- is.na(rho)
  base <- ifelse(
    is_na_rho,
    sprintf("n = %d (insufficient)", n),
    sprintf("rho = %.2f, p = %.2g\nn = %d", rho, p_value, n)
  )
  if (!is.null(slope)) {
    base <- ifelse(
      is_na_rho | is.na(slope),
      base,
      sprintf("rho = %.2f, p = %.2g\nslope = %.3f\nn = %d",
              rho, p_value, slope, n)
    )
  }
  base
}


# =============================================================================
# 1. Ortholog concordance — Q1
# =============================================================================

#' Scatter of human X vs mouse X, one point per ortholog gene pair, one panel
#' per metric. Identity line for reference. Per-panel Spearman summary.
#'
#' @param df        all_species dataframe.
#' @param metrics   Character vector of metric columns.
#'                  Default: all probing columns via fg_columns(df, "probing").
#' @param min_n     Minimum complete pairs per panel.
#' @param formatter Display formatter for facet labels (default format_col_name).
#' @return list(plot, table). Table: metric, metric_display, n, rho, p_value.
#' @export
ortholog_concordance_plot <- function(df,
                                      metrics   = fg_columns(df, "probing"),
                                      min_n     = 10,
                                      formatter = format_col_name) {
  paired <- pair_orthologs(df, metrics)
  paired_complete <- dplyr::filter(paired, !is.na(x_human), !is.na(x_mouse))
  
  summary_tbl <- spearman_summary(paired_complete, "x_human", "x_mouse",
                                  by = "metric", min_n = min_n) |>
    dplyr::mutate(metric_display = formatter(metric),
                  label = annotation_label(rho, p_value, n))
  
  # Order facets by |rho|, NAs last.
  metric_order <- summary_tbl |>
    dplyr::arrange(is.na(rho), dplyr::desc(abs(rho))) |>
    dplyr::pull(metric)
  paired_complete$metric <- factor(paired_complete$metric, levels = metric_order)
  summary_tbl$metric    <- factor(summary_tbl$metric,    levels = metric_order)
  
  annot <- paired_complete |>
    dplyr::group_by(metric) |>
    dplyr::summarise(
      x_pos = min(x_human, na.rm = TRUE),
      y_pos = max(x_mouse, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(summary_tbl, by = "metric")
  
  metric_labeller <- function(values) {
    formatter(as.character(values))
  }
  
  p <- ggplot2::ggplot(paired_complete,
                       ggplot2::aes(x = x_human, y = x_mouse)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.35, size = 1.2, colour = "#444444") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = unname(FEATURE_GROUP_COLOURS[["probing"]]),
                         linewidth = 0.7) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      hjust = 0, vjust = 1, size = 2.9, colour = "grey20",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free",
                        labeller = ggplot2::labeller(metric = metric_labeller)) +
    ggplot2::labs(
      title    = "Ortholog probing concordance: human vs mouse",
      subtitle = "One point per ortholog gene pair; dashed = identity",
      x        = "Human",
      y        = "Mouse"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey30", size = 10),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  table_out <- summary_tbl |>
    dplyr::transmute(metric = as.character(metric),
                     metric_display, n, rho, p_value)
  
  list(plot = p, table = tibble::as_tibble(table_out))
}


# =============================================================================
# 2. Per-species baseline — Q2
# =============================================================================

#' Faceted scatter of probing metric vs halflife, one panel per
#' (species, metric). Per-panel Spearman summary.
#'
#' @param df        all_species dataframe.
#' @param metrics   Character vector of metric columns.
#' @param regions   Character vector of region suffix tokens used to subset
#'                  `metrics` by their trailing `_<region>` token (e.g. "mrna",
#'                  or c("5utr","cds","3utr")). NULL (default) keeps all metrics.
#'                  Lets the caller split the panel grid by region.
#' @param min_n     Minimum n per panel.
#' @param formatter Display formatter for facet labels.
#' @return list(plot, table). Table: species, metric, metric_display,
#'         n, rho, p_value.
#' @export
probing_halflife_species_plot <- function(df,
                                          metrics   = fg_columns(df, "probing"),
                                          regions   = NULL,
                                          min_n     = 10,
                                          formatter = format_col_name) {
  # R5: guards.
  if (!"halflife" %in% names(df)) stop("halflife column missing from df")
  if (!"species"  %in% names(df)) stop("species column missing from df")
  
  # Optional region subset: keep metrics whose trailing _<region> token is in
  # `regions`. Suffix-based so it is robust to new compartments.
  if (!is.null(regions)) {
    region_pat <- paste0("[ _](", paste(regions, collapse = "|"), ")$")
    metrics <- metrics[grepl(region_pat, metrics)]
    if (length(metrics) == 0) {
      stop("no metrics match the requested regions: ",
           paste(regions, collapse = ", "))
    }
  }
  
  missing_metrics <- setdiff(metrics, names(df))
  if (length(missing_metrics)) {
    stop("metrics not found in df: ",
         paste(missing_metrics, collapse = ", "))
  }
  
  long <- df |>
    dplyr::select(species, halflife, dplyr::all_of(metrics)) |>
    tidyr::pivot_longer(dplyr::all_of(metrics),
                        names_to = "metric", values_to = "value") |>
    dplyr::filter(!is.na(value), !is.na(halflife))
  
  summary_tbl <- spearman_summary(long, "value", "halflife",
                                  by = c("species", "metric"),
                                  min_n = min_n) |>
    dplyr::mutate(metric_display = formatter(metric),
                  label = annotation_label(rho, p_value, n))
  
  # Order metric facets by absolute human rho (then mouse rho as tiebreaker).
  human_order <- summary_tbl |>
    dplyr::filter(species == "human") |>
    dplyr::arrange(is.na(rho), dplyr::desc(abs(rho))) |>
    dplyr::pull(metric)
  metric_order <- c(human_order,
                    setdiff(unique(summary_tbl$metric), human_order))
  long$metric        <- factor(long$metric,        levels = metric_order)
  summary_tbl$metric <- factor(summary_tbl$metric, levels = metric_order)
  
  annot <- long |>
    dplyr::group_by(species, metric) |>
    dplyr::summarise(
      x_pos = min(value, na.rm = TRUE),
      y_pos = max(halflife, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(summary_tbl, by = c("species", "metric"))
  
  metric_labeller <- function(values) {
    formatter(as.character(values))
  }
  
  species_colours <- c(human = "#1F78B4", mouse = "#E66101")
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = value, y = halflife,
                                          colour = species)) +
    ggplot2::geom_point(alpha = 0.2, size = 0.9) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         linewidth = 0.8) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      hjust = 0, vjust = 1, size = 2.7, colour = "grey20",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_grid(species ~ metric, scales = "free",
                        labeller = ggplot2::labeller(metric = metric_labeller)) +
    ggplot2::scale_colour_manual(values = species_colours, name = "Species") +
    ggplot2::labs(
      title    = "Probing metric vs half-life, per species",
      subtitle = "Spearman rho computed within each species; metrics ordered by |rho| in human",
      x        = "Probing metric value",
      y        = format_col_name("halflife")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey30", size = 10),
      legend.position  = "none",  # species shown on row strip
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  table_out <- summary_tbl |>
    dplyr::transmute(species, metric = as.character(metric),
                     metric_display, n, rho, p_value)
  
  list(plot = p, table = tibble::as_tibble(table_out))
}


# =============================================================================
# 3. Ortholog Δ-Δ — Q3
# =============================================================================

#' Within-ortholog Δ-Δ scatter: Δprobing (human − mouse) vs Δhalflife
#' (human − mouse), faceted by metric. Per-panel Spearman + linear-model slope.
#'
#' Differencing on ortholog pairs controls for gene-level confounders
#' (function, expression class, length, GC), so this is the strongest test
#' that a probing → halflife relationship is mechanistically conserved.
#'
#' @param df        all_species dataframe.
#' @param metrics   Default fg_columns(df, "probing").
#' @param min_n     Minimum ortholog pairs per panel.
#' @param formatter Display formatter for facet labels.
#' @return list(plot, table). Table: metric, metric_display, n, rho, p_value,
#'         lm_slope, lm_intercept, lm_slope_p.
#' @export
ortholog_delta_delta_plot <- function(df,
                                      metrics   = fg_columns(df, "probing"),
                                      min_n     = 10,
                                      formatter = format_col_name) {
  paired <- pair_orthologs(df, metrics)
  paired_complete <- dplyr::filter(paired, !is.na(dx), !is.na(dhl))
  
  spear_tbl <- spearman_summary(paired_complete, "dx", "dhl",
                                by = "metric", min_n = min_n)
  
  lm_tbl <- paired_complete |>
    dplyr::group_by(metric) |>
    dplyr::summarise(
      .lm = list({
        if (dplyr::n() < min_n) {
          list(slope = NA_real_, intercept = NA_real_, slope_p = NA_real_)
        } else {
          fit <- stats::lm(dhl ~ dx)
          co  <- summary(fit)$coefficients
          list(
            slope     = unname(co[2, "Estimate"]),
            intercept = unname(co[1, "Estimate"]),
            slope_p   = unname(co[2, "Pr(>|t|)"])
          )
        }
      }),
      .groups = "drop"
    ) |>
    tidyr::unnest_wider(.lm) |>
    dplyr::rename(lm_slope = slope, lm_intercept = intercept,
                  lm_slope_p = slope_p)
  
  summary_tbl <- spear_tbl |>
    dplyr::left_join(lm_tbl, by = "metric") |>
    dplyr::mutate(metric_display = formatter(metric),
                  label = annotation_label(rho, p_value, n, slope = lm_slope))
  
  metric_order <- summary_tbl |>
    dplyr::arrange(is.na(rho), dplyr::desc(abs(rho))) |>
    dplyr::pull(metric)
  paired_complete$metric <- factor(paired_complete$metric, levels = metric_order)
  summary_tbl$metric    <- factor(summary_tbl$metric,    levels = metric_order)
  
  annot <- paired_complete |>
    dplyr::group_by(metric) |>
    dplyr::summarise(
      x_pos = min(dx, na.rm = TRUE),
      y_pos = max(dhl, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(summary_tbl, by = "metric")
  
  metric_labeller <- function(values) {
    formatter(as.character(values))
  }
  
  probing_col <- unname(FEATURE_GROUP_COLOURS[["probing"]])
  
  p <- ggplot2::ggplot(paired_complete, ggplot2::aes(x = dx, y = dhl)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.4, size = 1.2, colour = "#444444") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         colour = probing_col, fill = probing_col,
                         alpha = 0.15, linewidth = 0.7) +
    ggplot2::geom_text(
      data = annot,
      ggplot2::aes(x = x_pos, y = y_pos, label = label),
      hjust = 0, vjust = 1, size = 2.8, colour = "grey20",
      inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~ metric, scales = "free",
                        labeller = ggplot2::labeller(metric = metric_labeller)) +
    ggplot2::labs(
      title    = "Ortholog Delta-Delta: does Delta probing track Delta half-life?",
      subtitle = "Delta = human - mouse; one point per ortholog gene pair",
      x        = "Delta probing metric (human - mouse)",
      y        = "Delta half-life (human - mouse)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey30", size = 10),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text       = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  table_out <- summary_tbl |>
    dplyr::transmute(metric = as.character(metric), metric_display,
                     n, rho, p_value, lm_slope, lm_intercept, lm_slope_p)
  
  list(plot = p, table = tibble::as_tibble(table_out))
}


# =============================================================================
# Runner block — executes only when this file is sourced top-to-bottom.
# =============================================================================

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  
  all_species <- build_all()
  
  out_dir <- file.path(OUTPUT_DIR, "plots", "cross_species")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  q1 <- ortholog_concordance_plot(all_species)
  q2_mrna    <- probing_halflife_species_plot(all_species, regions = "mrna")
  q2_regions <- probing_halflife_species_plot(all_species,
                                              regions = c("5utr", "cds", "3utr"))
  q3 <- ortholog_delta_delta_plot(all_species)
  
  ggplot2::ggsave(file.path(out_dir, "ortholog_concordance.pdf"),
                  q1$plot, width = 10, height = 7)
  ggplot2::ggsave(file.path(out_dir, "probing_halflife_species_mrna.pdf"),
                  q2_mrna$plot, width = 7, height = 6)
  ggplot2::ggsave(file.path(out_dir, "probing_halflife_species_regions.pdf"),
                  q2_regions$plot, width = 14, height = 6)
  ggplot2::ggsave(file.path(out_dir, "ortholog_delta_delta.pdf"),
                  q3$plot, width = 10, height = 7)
  
  readr::write_tsv(q1$table,
                   file.path(out_dir, "ortholog_concordance.tsv"))
  readr::write_tsv(q2_mrna$table,
                   file.path(out_dir, "probing_halflife_species_mrna.tsv"))
  readr::write_tsv(q2_regions$table,
                   file.path(out_dir, "probing_halflife_species_regions.tsv"))
  readr::write_tsv(q3$table,
                   file.path(out_dir, "ortholog_delta_delta.tsv"))
  
  message("Wrote 4 plots + 4 tables to: ", out_dir)
}