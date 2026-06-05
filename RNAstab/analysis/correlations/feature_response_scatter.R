# =============================================================================
# Feature × two-response correlation scatter
# =============================================================================
# Generalised refactor of the legacy corr_te_halflife_plot. Each point is a
# feature, positioned by its correlation with response_x against response_y.
#
# Quadrant interpretation:
#   top-right / bottom-left  — concordant predictors (same sign for both)
#   top-left / bottom-right  — discordant (e.g. translation-vs-stability tradeoffs)
#   near origin              — uninformative for both
#   far along one axis only  — response-specific predictors
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_response_scatter.R")
#   df  <- build_dataset("human")
#
#   # Default: TE vs halflife, all features, per-column granularity
#   out <- feature_response_scatter(df)
#   print(out$plot)
#
#   # Cleaner: collapse to one point per (group, region), filter noise:
#   out <- feature_response_scatter(
#     df, collapse = "region", noise_filter = 0.1
#   )
#
#   # Restrict to structure features only:
#   out <- feature_response_scatter(df, groups = "structure")
#
#   # Saluki diagnostic: which features explain Saluki's residuals?
#   out <- feature_response_scatter(
#     df,
#     response_x = "saluki_prediction",
#     response_y = "prediction_difference",
#     exclude    = c("^halflife$")
#   )
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(tibble)
})


#' Scatter of feature correlations against two responses.
#'
#' @param df             Dataframe from build_dataset() or build_all().
#' @param response_x     Character. Column for the x-axis correlation.
#' @param response_y     Character. Column for the y-axis correlation.
#' @param method         Correlation method (default "spearman").
#' @param groups         Character vector of FEATURE_PATTERNS keys, SUPERGROUPS
#'                       names, and/or GROUP_BUNDLES names. NULL (default) =
#'                       all groups. "structure" expands to all structure
#'                       groups; a bundle expands to its groups + pick/drop.
#' @param pick           Named list: group key -> columns to keep (allow-list,
#'                       caller order). New columns added to the group later
#'                       stay out until named. Merges with any bundle pick
#'                       (caller wins per group key).
#' @param drop           Named list: group key -> columns to remove from the
#'                       otherwise-whole group. New columns added later are
#'                       included. Merges with any bundle drop (caller wins).
#' @param collapse       "none" (one point per column, default), "region" (one
#'                       point per group × region — median r), or "group" (one
#'                       point per group — median r across all members).
#' @param noise_filter   Numeric. Drop points with distance-from-origin below
#'                       this. 0 (default) = no filter. Try 0.1 to declutter.
#' @param label_quantile Numeric in [0, 1]. Label the top (1 - q) fraction by
#'                       distance from origin. Default 0.9 = label top 10%.
#' @param standalones    Character vector of standalone columns to include
#'                       even though they're not in any group. Default includes
#'                       cai, translation_efficiency, expression.
#' @param exclude        Regex patterns to exclude. Default = R10 derived-
#'                       prediction set, unless response_x or response_y is
#'                       one of them (in which case it's auto-removed from
#'                       the exclude list).
#' @param min_n          Minimum non-NA pairs to compute a correlation.
#' @param formatter      Display formatter (default format_col_name).
#' @param palette        Named colour vector keyed by group/standalone name.
#'                       Default FEATURE_GROUP_COLOURS.
#' @param shapes         Named shape vector keyed by region token.
#'                       Default REGION_SHAPES.
#' @return list(plot, table). Table columns:
#'   species, variable, group, supergroup, region, n_x, correlation_x,
#'   p_value_x, q_x, n_y, correlation_y, p_value_y, q_y,
#'   distance_from_origin, labelled.
#' @export
feature_response_scatter <- function(df,
                                     response_x     = "translation_efficiency",
                                     response_y     = "halflife",
                                     method         = c("spearman", "pearson",
                                                        "kendall"),
                                     groups         = NULL,
                                     pick           = list(),
                                     drop           = list(),
                                     collapse       = c("none", "region",
                                                        "group"),
                                     noise_filter   = 0,
                                     label_quantile = 0.9,
                                     standalones    = c("cai",
                                                        "translation_efficiency",
                                                        "expression"),
                                     exclude        = c("^saluki_prediction$",
                                                        "^prediction_difference$"),
                                     min_n          = 30,
                                     formatter      = format_col_name,
                                     palette        = NULL,
                                     shapes         = NULL) {

  method   <- match.arg(method)
  collapse <- match.arg(collapse)
  if (is.null(palette)) palette <- FEATURE_GROUP_COLOURS
  if (is.null(shapes))  shapes  <- REGION_SHAPES

  # --- R5: guard ----------------------------------------------------------
  if (!response_x %in% names(df)) stop("response_x '", response_x, "' not in df")
  if (!response_y %in% names(df)) stop("response_y '", response_y, "' not in df")
  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }
  if (response_x == response_y) {
    stop("response_x and response_y must differ")
  }
  if (label_quantile < 0 || label_quantile >= 1) {
    stop("label_quantile must be in [0, 1)")
  }

  # Don't exclude a response if the user is correlating against it
  exclude <- exclude[!vapply(exclude, function(rgx) {
    grepl(rgx, response_x) || grepl(rgx, response_y)
  }, logical(1))]

  # --- Enumerate candidate columns ----------------------------------------
  sel      <- resolve_selection(groups, pick, drop)
  expanded <- sel$groups

  # Build a column -> group attribution map. First-match wins, so aliases
  # (already removed from FEATURE_PATTERNS) don't cause double-assignment.
  col_to_group <- list()
  for (g in expanded) {
    # Refine via the shared helper so bundle- and caller-supplied pick/drop
    # (already merged by resolve_selection) apply identically to select_features.
    cols <- refine_group_columns(fg_columns(df, g), sel$pick[[g]], sel$drop[[g]])
    for (co in cols) {
      if (is.null(col_to_group[[co]])) col_to_group[[co]] <- g
    }
  }

  # Add standalone columns of interest (not already in a group)
  for (co in standalones) {
    if (co %in% names(df) && is.null(col_to_group[[co]])) {
      col_to_group[[co]] <- co        # standalone: group key = column name
    }
  }

  candidates <- names(col_to_group)

  # Apply exclusions
  for (rgx in exclude) {
    candidates <- candidates[!vapply(candidates, function(c) grepl(rgx, c),
                                     logical(1))]
  }

  # Never include the responses themselves
  candidates <- setdiff(candidates, c(response_x, response_y))

  if (length(candidates) == 0) {
    stop("No candidate features after filtering — check `groups` and `exclude`")
  }

  # --- Per-species correlation computation --------------------------------
  has_species <- length(unique(df$species)) > 1

  compute_one <- function(sub, sp_label) {
    purrr::map_dfr(candidates, function(co) {
      v  <- sub[[co]]
      vx <- sub[[response_x]]
      vy <- sub[[response_y]]

      ok_x <- !is.na(v) & !is.na(vx)
      ok_y <- !is.na(v) & !is.na(vy)

      if (sum(ok_x) < min_n || sum(ok_y) < min_n) return(tibble::tibble())
      if (length(unique(v[ok_x])) < 2 ||
          length(unique(v[ok_y])) < 2) return(tibble::tibble())

      ct_x <- suppressWarnings(stats::cor.test(
        v[ok_x], vx[ok_x], method = method, exact = FALSE
      ))
      ct_y <- suppressWarnings(stats::cor.test(
        v[ok_y], vy[ok_y], method = method, exact = FALSE
      ))

      tibble::tibble(
        species       = sp_label,
        variable      = co,
        group         = col_to_group[[co]],
        n_x           = sum(ok_x),
        correlation_x = unname(ct_x$estimate),
        p_value_x     = ct_x$p.value,
        n_y           = sum(ok_y),
        correlation_y = unname(ct_y$estimate),
        p_value_y     = ct_y$p.value
      )
    })
  }

  result <- if (has_species) {
    purrr::map_dfr(unique(df$species), function(sp) {
      compute_one(df |> dplyr::filter(species == sp), sp)
    })
  } else {
    compute_one(df, unique(df$species)[1])
  }

  if (nrow(result) == 0) {
    stop("No correlations computed — try lowering min_n or check coverage")
  }

  # --- Region extraction (REGIONS-aware) ----------------------------------
  result <- result |>
    dplyr::mutate(
      region = vapply(variable, function(co) {
        tokens <- strsplit(co, "_", fixed = TRUE)[[1]]
        last <- tokens[length(tokens)]
        if (last %in% REGIONS) last else "none"
      }, character(1)),
      supergroup = dplyr::coalesce(supergroup_of(group), "other")
    )

  # --- BH q-values, per species, per response axis ------------------------
  result <- result |>
    dplyr::group_by(species) |>
    dplyr::mutate(
      q_x = stats::p.adjust(p_value_x, method = "BH"),
      q_y = stats::p.adjust(p_value_y, method = "BH")
    ) |>
    dplyr::ungroup()

  # --- Collapse if requested ----------------------------------------------
  if (collapse != "none") {
    grp_cols <- switch(collapse,
                       region = c("species", "group", "region"),
                       group  = c("species", "group"))

    result <- result |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) |>
      dplyr::summarise(
        n_variables       = dplyr::n(),
        representative    = variable[which.max(abs(correlation_x) +
                                               abs(correlation_y))],
        correlation_x     = stats::median(correlation_x, na.rm = TRUE),
        correlation_y     = stats::median(correlation_y, na.rm = TRUE),
        n_x               = stats::median(n_x, na.rm = TRUE),
        n_y               = stats::median(n_y, na.rm = TRUE),
        q_x               = stats::median(q_x, na.rm = TRUE),
        q_y               = stats::median(q_y, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        # Variable name: synthetic, useful for the table identifier
        variable   = if (collapse == "region") {
          paste(group, region, sep = "__")
        } else {
          group
        },
        # Recompute supergroup (collapse may have dropped it)
        supergroup = dplyr::coalesce(supergroup_of(group), "other"),
        # Drop region in group-only collapse
        region     = if (collapse == "group") "none" else region
      )
  }

  # --- Distance + noise filter --------------------------------------------
  result <- result |>
    dplyr::mutate(
      distance_from_origin = sqrt(correlation_x^2 + correlation_y^2)
    ) |>
    dplyr::filter(distance_from_origin >= noise_filter)

  if (nrow(result) == 0) {
    stop("All points filtered out — noise_filter too high?")
  }

  # --- Label selection by quantile ----------------------------------------
  # Per species: threshold at the requested quantile of distance.
  result <- result |>
    dplyr::group_by(species) |>
    dplyr::mutate(
      .thr     = stats::quantile(distance_from_origin,
                                 probs = label_quantile, na.rm = TRUE),
      labelled = distance_from_origin >= .thr
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.thr)

  # --- Display labels (R4) -------------------------------------------------
  result <- result |>
    dplyr::mutate(
      display_label = dplyr::case_when(
        collapse == "group"  ~ format_group_name(group),
        collapse == "region" ~ paste0(
          format_group_name(group),
          ifelse(region == "none", "",
                 paste0(" — ", format_col_name(region)))
        ),
        TRUE                 ~ formatter(variable)
      ),
      label_text = ifelse(labelled, display_label, NA_character_)
    )

  # --- Build the plot ------------------------------------------------------
  # Use group as the colour aesthetic; legend label via format_group_name
  # (or format_col_name for standalones). Region is shape.
  legend_labeller <- function(g) {
    ifelse(g %in% names(FEATURE_PATTERNS) |
             g %in% c("structure", "intrinsic", "splicing",
                      "regulatory", "other"),
           format_group_name(g),
           formatter(g))
  }

  # Make the categorical axes factors so legend ordering is consistent
  group_order <- intersect(names(palette), unique(result$group))
  result$group_f <- factor(result$group, levels = group_order)
  shape_order <- intersect(names(shapes), unique(result$region))
  result$region_f <- factor(result$region, levels = shape_order)

  axis_x_lab <- sprintf("%s correlation with %s",
                        tools::toTitleCase(method), formatter(response_x))
  axis_y_lab <- sprintf("%s correlation with %s",
                        tools::toTitleCase(method), formatter(response_y))

  title <- sprintf("Feature correlations: %s vs %s",
                   formatter(response_x), formatter(response_y))
  subtitle_bits <- character()
  if (collapse != "none") {
    subtitle_bits <- c(subtitle_bits,
                       paste0("collapsed by ", collapse,
                              " (median ", method, ")"))
  } else {
    subtitle_bits <- c(subtitle_bits, paste0("per-column ", method))
  }
  if (noise_filter > 0) {
    subtitle_bits <- c(subtitle_bits,
                       sprintf("noise filter |r| >= %.2f", noise_filter))
  }
  subtitle_bits <- c(subtitle_bits,
                     sprintf("top %d%% labelled",
                             round((1 - label_quantile) * 100)))

  p <- ggplot2::ggplot(
    result,
    ggplot2::aes(x = correlation_x, y = correlation_y)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_point(
      ggplot2::aes(colour = group_f, shape = region_f),
      size = 3, alpha = 0.85, stroke = 0.7
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label_text, colour = group_f),
      size = 3, max.overlaps = 50,
      box.padding = 0.5, min.segment.length = 0,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = palette,
      labels = legend_labeller,
      name   = "Feature group",
      drop   = TRUE
    ) +
    ggplot2::scale_shape_manual(
      values = shapes,
      labels = function(r) ifelse(r == "none", "transcript-level",
                                  formatter(r)),
      name   = "Region",
      drop   = TRUE
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = paste(subtitle_bits, collapse = " · "),
      x        = axis_x_lab,
      y        = axis_y_lab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30"),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(face = "bold"),
      legend.text      = ggplot2::element_text(size = 8),
      legend.key.size  = ggplot2::unit(0.8, "lines"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1, ncol = 1,
                                     override.aes = list(size = 3, alpha = 1)),
      shape  = ggplot2::guide_legend(order = 2,
                                     override.aes = list(size = 3))
    )

  if (has_species) {
    p <- p + ggplot2::facet_wrap(~ species)
  }

  # --- R9: return table without the plot-only helper columns -------------
  table_out <- result |>
    dplyr::select(dplyr::any_of(c(
      "species", "variable", "group", "supergroup", "region",
      "n_variables", "representative",
      "n_x", "correlation_x", "p_value_x", "q_x",
      "n_y", "correlation_y", "p_value_y", "q_y",
      "distance_from_origin", "labelled", "display_label"
    )))

  list(plot = p, table = table_out)
}


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  df <- build_dataset("human")

  # --- R8: outputs go under OUTPUT_DIR ------------------------------------
  dir.create(file.path(OUTPUT_DIR, "plots"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"),
             showWarnings = FALSE, recursive = TRUE)

  # Default view: TE vs halflife, all features, per-column
  out_full <- feature_response_scatter(df)
  print(out_full$plot)
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "plots",
              "feature_response_scatter_te_vs_halflife_full.jpg"),
    plot = out_full$plot,
    width = 280, height = 200, units = "mm", dpi = 300
  )
  write.csv(
    out_full$table,
    file.path(OUTPUT_DIR, "tables",
              "feature_response_scatter_te_vs_halflife_full.csv"),
    row.names = FALSE
  )

  # Cleaner view: collapse to (group, region), filter noise
  out_clean <- feature_response_scatter(
    df,
    collapse     = "region",
    noise_filter = 0.05
  )
  print(out_clean$plot)
  ggplot2::ggsave(
    file.path(OUTPUT_DIR, "plots",
              "feature_response_scatter_te_vs_halflife_collapsed.jpg"),
    plot = out_clean$plot,
    width = 250, height = 180, units = "mm", dpi = 300
  )
  write.csv(
    out_clean$table,
    file.path(OUTPUT_DIR, "tables",
              "feature_response_scatter_te_vs_halflife_collapsed.csv"),
    row.names = FALSE
  )

  message("\nFeature response scatter complete:")
  message("  ", file.path(OUTPUT_DIR, "plots"),
          "/feature_response_scatter_*.jpg")
  message("  ", file.path(OUTPUT_DIR, "tables"),
          "/feature_response_scatter_*.csv")
}
