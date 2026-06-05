# =============================================================================
# Feature correlation dotplot with confidence intervals
# =============================================================================
# Generalised refactor of the legacy plotting_corr_w_CIs script.
#
# Each x-axis tick is a metric "stem" (region stripped); each dodged dot is
# the correlation of one region's version of that metric with a chosen
# response, with confidence-interval error bars. Faceted by SUPERGROUPS so
# related categories sit together. The CI is essential context the scatter
# plot in feature_response_scatter doesn't show — same correlation with
# different sample sizes carries different weight.
#
# Generalised vs the original:
#   * Any response, not just halflife.
#   * Multi-species ready (faceted as species x supergroup).
#   * Drives off FEATURE_PATTERNS + SUPERGROUPS for categorisation.
#   * REGIONS-aware extraction (handles last100, utrpair, etc. correctly).
#   * Optional |r| transformation (default ON, matching the original).
#   * Adjustable significance threshold line.
#
# --- v4 region vocabulary ---------------------------------------------------
# Every splicing and regulatory column carries a real region token as its
# LAST token:
#   * junctions   -> junctions_count_<5utr|cds|3utr>
#   * distances   -> eej_dist_<first|last>_<start|stop>   (start/stop are
#                                                          real regions)
#   * architecture-> *_mrna                               (whole-transcript)
#   * uorfs       -> uorf_*_mrna                           (whole-transcript)
#   * nmd         -> nmd_*_mrna                            (single model)
# so the whole splicing supergroup and the transcript-level uORF/architecture
# columns are picked up by the region-token filter automatically. Two further
# tiers are handled explicitly:
#   * Tier 1 — reserved single-token scalars (cai, expression,
#     translation_efficiency, orfexondensity): genuinely region-less. Mapped
#     to the `mrna` region via the `standalones` argument so they appear in
#     the "other" supergroup facet rather than being dropped.
#   * Any column that still has no region token is reported by the diagnostic
#     block and skipped (it cannot sit on the region-dodged axis).
#
# CI method:
#   * Pearson  : cor.test() built-in conf.int (exact for normal data).
#   * Spearman : Fisher's Z back-transform with Bonett-Wright SE
#                = 1.06 / sqrt(n - 3). Accurate for n > ~30.
#   * Kendall  : Fisher's Z back-transform with SE = sqrt(0.437 / (n - 4)).
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/correlations/feature_correlation_dotplot.R")
#   df  <- build_dataset("human")
#
#   # Default: halflife
#   out <- feature_correlation_dotplot(df)
#   print(out$plot)
#
#   # Translation efficiency
#   out <- feature_correlation_dotplot(df, response = "translation_efficiency")
#
#   # Signed correlations instead of absolute
#   out <- feature_correlation_dotplot(df, absolute = FALSE)
#
#   # Restrict to structure features, drop low-|r| noise
#   out <- feature_correlation_dotplot(df, groups = "structure",
#                                      min_abs_correlation = 0.05)
#
#   # Drop the transcript-level standalone scalars
#   out <- feature_correlation_dotplot(df, standalones = character())
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


# -----------------------------------------------------------------------------
# Helper: correlation with CI for one (x, y) pair
# -----------------------------------------------------------------------------

#' Compute a correlation and a Fisher's-Z confidence interval.
#' Returns NA estimate/CI if too few non-NA pairs.
#'
#' @param x,y      Numeric vectors of equal length.
#' @param method   "spearman", "pearson", or "kendall".
#' @param conf     Confidence level (default 0.95).
#' @param min_n    Minimum non-NA pairs to attempt.
#' @return list(estimate, conf.low, conf.high, p.value, n).
correlation_with_ci <- function(x, y,
                                method = c("spearman", "pearson", "kendall"),
                                conf   = 0.95,
                                min_n  = 30) {
  method <- match.arg(method)
  ok <- !is.na(x) & !is.na(y)
  n  <- sum(ok)

  if (n < min_n ||
      length(unique(x[ok])) < 2 ||
      length(unique(y[ok])) < 2) {
    return(list(estimate = NA_real_, conf.low = NA_real_,
                conf.high = NA_real_, p.value = NA_real_, n = n))
  }

  ct <- suppressWarnings(
    stats::cor.test(x[ok], y[ok], method = method, exact = FALSE,
                    conf.level = conf)
  )

  r <- unname(ct$estimate)

  # Pearson's cor.test gives a conf.int natively. For Spearman / Kendall we
  # build one via Fisher's Z with the appropriate variance correction.
  if (method == "pearson" && !is.null(ct$conf.int)) {
    return(list(
      estimate  = r,
      conf.low  = ct$conf.int[1],
      conf.high = ct$conf.int[2],
      p.value   = ct$p.value,
      n         = n
    ))
  }

  z   <- atanh(r)
  se  <- switch(
    method,
    spearman = 1.06 / sqrt(n - 3),
    kendall  = sqrt(0.437 / (n - 4))
  )
  zc  <- stats::qnorm(1 - (1 - conf) / 2)
  list(
    estimate  = r,
    conf.low  = tanh(z - zc * se),
    conf.high = tanh(z + zc * se),
    p.value   = ct$p.value,
    n         = n
  )
}


# -----------------------------------------------------------------------------
# Main plot function
# -----------------------------------------------------------------------------

#' Feature correlations with a chosen response, with confidence intervals,
#' faceted by supergroup, regions as dodged points.
#'
#' @param df                 Dataframe from build_dataset() or build_all().
#' @param response           Character. Response column (default "halflife").
#' @param method             Correlation method (default "spearman").
#' @param groups             Character vector of FEATURE_PATTERNS keys and/or
#'                           SUPERGROUPS names. NULL (default) = all groups.
#' @param pick               Named list: group key -> columns to keep (allow-
#'                           list; caller order honoured for the within-group
#'                           sequence, though the plot re-orders by max |r|).
#'                           New columns added to the group later stay out
#'                           until named. Use for a small fixed subset of an
#'                           open-ended family.
#' @param drop               Named list: group key -> columns to remove from
#'                           the otherwise-whole group. New columns added later
#'                           are included. Use for "the family minus a couple
#'                           of noisy members".
#' @param standalones        Character vector of reserved single-token scalar
#'                           columns to include. These have no region suffix
#'                           and are mapped to the `mrna` region so they
#'                           appear on the dotplot (in the "other" supergroup
#'                           facet). Default = cai,
#'                           translation_efficiency, orfexondensity. Pass
#'                           character() to exclude them.
#' @param absolute           Logical. If TRUE (default) plot |correlation|
#'                           and transform the CI accordingly; if FALSE plot
#'                           signed correlation.
#' @param min_abs_correlation Numeric. Drop points where |r| is below this.
#'                           0 (default) = no filter.
#' @param sig_threshold      Numeric or NULL. Reference horizontal line at
#'                           this y value (e.g. 0.02 matching the original).
#'                           NULL = no line.
#' @param conf               Confidence level for the CI (default 0.95).
#' @param min_n              Minimum non-NA pairs to compute a correlation.
#' @param formatter          Display formatter for the response label
#'                           (default format_col_name).
#' @param region_colours     Named colour vector. Default REGION_COLOURS.
#' @param region_shapes      Named shape vector. Default REGION_SHAPES.
#' @param top_n_per_group  Named list. For each named FEATURE_PATTERNS group,
#'                         restrict the plot to the top-N metric stems within
#'                         that group, ranked by max |correlation| across
#'                         regions and species. Use for high-cardinality
#'                         groups like codon_freqs (~64) and aa_freqs (~20),
#'                         which are excluded from the default expansion and
#'                         must be passed explicitly via `groups`.
#'                         Example: top_n_per_group = list(codon_freqs = 10,
#'                                                         aa_freqs    = 5)
#' @return list(plot, table). Table columns: species, variable, group,
#'   supergroup, metric_stem, metric_display, region, n, correlation,
#'   conf.low, conf.high, p_value, q_value, correlation_abs, conf.low_abs,
#'   conf.high_abs.
#' @export
feature_correlation_dotplot <- function(df,
                                        response             = "halflife",
                                        method               = c("spearman",
                                                                 "pearson",
                                                                 "kendall"),
                                        groups               = NULL,
                                        pick                 = list(),
                                        drop                 = list(),
                                        standalones          = c("cai",
                                                                 "translation_efficiency",
                                                                 "orfexondensity"),
                                        absolute             = TRUE,
                                        min_abs_correlation  = 0,
                                        sig_threshold        = NULL,
                                        conf                 = 0.95,
                                        min_n                = 30,
                                        formatter            = format_col_name,
                                        region_colours       = NULL,
                                        region_shapes        = NULL,
                                        top_n_per_group = list()) {

  method <- match.arg(method)
  if (is.null(region_colours)) region_colours <- REGION_COLOURS
  if (is.null(region_shapes))  region_shapes  <- REGION_SHAPES

  # --- R5: guard ----------------------------------------------------------
  if (!response  %in% names(df)) stop("response '", response, "' not in df")
  if (!"species" %in% names(df)) {
    stop("species column missing - pipeline invariant violated")
  }

  # --- Enumerate region-bearing columns -----------------------------------
  # v4: every splicing/regulatory column ends in a real region token (mrna
  # for whole-transcript scalars, start/stop for EEJ distances, 5utr/cds/3utr
  # for junction counts), so they are picked up here automatically. Only
  # genuinely malformed columns (no region token at all) fall through to
  # `dropped`.
  sel      <- resolve_selection(groups, pick, drop)
  expanded <- sel$groups

  col_to_group <- list()   # column -> FEATURE_PATTERNS key (or standalone name)
  col_region   <- list()   # column -> region token (real or pseudo)
  col_stem     <- list()   # column -> metric stem (column minus region token)
  dropped      <- character()

  for (g in expanded) {
    # Refine via the shared helper so bundle- and caller-supplied pick/drop
    # (already merged by resolve_selection) apply identically to select_features.
    cols <- refine_group_columns(fg_columns(df, g), sel$pick[[g]], sel$drop[[g]])

    for (co in cols) {
      if (!is.null(col_to_group[[co]])) next
      tokens <- strsplit(co, "_", fixed = TRUE)[[1]]
      last   <- tokens[length(tokens)]
      if (last %in% REGIONS && length(tokens) > 1) {
        col_to_group[[co]] <- g
        col_region[[co]]   <- last
        col_stem[[co]]     <- paste(tokens[-length(tokens)], collapse = "_")
      } else {
        # No region token — cannot sit on the region-dodged axis.
        dropped <- c(dropped, co)
      }
    }
  }

  # --- v4 tier 1: reserved single-token scalar columns --------------------
  # cai / expression / translation_efficiency / orfexondensity are genuinely
  # region-less. Map them to the `mrna` region so they appear on the dotplot
  # (in the "other" supergroup facet) rather than being silently excluded.
  # supergroup_of() returns NA for these -> coalesced to "other".
  for (co in standalones) {
    if (co %in% names(df) && is.null(col_to_group[[co]])) {
      col_to_group[[co]] <- co            # standalone: group key = column name
      col_region[[co]]   <- "mrna"
      col_stem[[co]]     <- co
    }
  }

  candidates <- setdiff(names(col_to_group), response)
  if (length(candidates) == 0) {
    stop("No plottable columns after filtering - check `groups` / `standalones`")
  }

  # --- Diagnostic: report columns dropped for lacking a region token ------
  # Turns the silent exclusion into a visible per-group report. After the v3
  # rename this should normally be empty; a non-empty report flags a column
  # whose loader still violates the region-suffix-last invariant.
  dropped <- setdiff(unique(dropped), names(col_to_group))
  if (length(dropped) > 0) {
    message("feature_correlation_dotplot: dropped ", length(dropped),
            " column(s) with no region token (not plottable here):")
    drop_by_group <- vapply(expanded, function(g) {
      sum(fg_columns(df, g) %in% dropped)
    }, integer(1))
    for (g in expanded[drop_by_group > 0]) {
      message("  ", g, ": ", drop_by_group[[which(expanded == g)]],
              " column(s)")
    }
    message("  -> use top_n_response_correlations() or ",
            "feature_response_scatter() for these")
  }

  # --- Per-species correlation + CI ---------------------------------------
  has_species <- length(unique(df$species)) > 1

  compute_one <- function(sub, sp_label) {
    purrr::map_dfr(candidates, function(co) {
      r <- correlation_with_ci(
        sub[[co]], sub[[response]],
        method = method, conf = conf, min_n = min_n
      )
      if (is.na(r$estimate)) return(tibble::tibble())

      tibble::tibble(
        species     = sp_label,
        variable    = co,
        group       = col_to_group[[co]],
        region      = col_region[[co]],
        metric_stem = col_stem[[co]],
        n           = r$n,
        correlation = r$estimate,
        conf.low    = r$conf.low,
        conf.high   = r$conf.high,
        p_value     = r$p.value
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
    stop("No correlations computed - try lowering min_n or check coverage")
  }

  # --- Attach supergroup + BH q -------------------------------------------
  result <- result |>
    dplyr::mutate(
      supergroup = dplyr::coalesce(supergroup_of(group), "other")
    ) |>
    dplyr::group_by(species) |>
    dplyr::mutate(q_value = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup()

  # --- |r| transform (preserving CI semantics) ----------------------------
  # When the CI straddles zero, |r| CI runs [0, max(|low|, |high|)].
  # When the CI is fully negative, |r| CI is [|high|, |low|].
  # When fully positive, |r| CI is unchanged.
  result <- result |>
    dplyr::mutate(
      correlation_abs = abs(correlation),
      conf.low_abs = dplyr::case_when(
        conf.low < 0 & conf.high > 0 ~ 0,
        conf.low < 0 & conf.high < 0 ~ abs(conf.high),
        TRUE                         ~ conf.low
      ),
      conf.high_abs = dplyr::case_when(
        conf.low < 0 & conf.high > 0 ~ pmax(abs(conf.low), abs(conf.high)),
        conf.low < 0 & conf.high < 0 ~ abs(conf.low),
        TRUE                         ~ conf.high
      )
    )

  # --- Top-N per group filter (for high-cardinality groups) ---------------
  # Rank by metric_stem so all regions of a stem stay together (no half-
  # plotted stems). Ranking aggregates across species too — keeps the
  # selected stems consistent between species panels.
  if (length(top_n_per_group) > 0) {
    for (g in names(top_n_per_group)) {
      n_keep <- top_n_per_group[[g]]
      if (!any(result$group == g)) next
      
      keep_stems <- result |>
        dplyr::filter(group == g) |>
        dplyr::group_by(metric_stem) |>
        dplyr::summarise(max_r = max(correlation_abs, na.rm = TRUE),
                         .groups = "drop") |>
        dplyr::arrange(dplyr::desc(max_r)) |>
        dplyr::slice_head(n = n_keep) |>
        dplyr::pull(metric_stem)
      
      result <- result |>
        dplyr::filter(group != g | metric_stem %in% keep_stems)
    }
  }
  
  # --- Filter by |r| threshold --------------------------------------------
  if (min_abs_correlation > 0) {
    result <- result |>
      dplyr::filter(correlation_abs >= min_abs_correlation)
    if (nrow(result) == 0) {
      stop("All points filtered out — min_abs_correlation too high?")
    }
  }

  # --- Pick the y aesthetic + CI columns based on absolute flag -----------
  if (absolute) {
    result$.y    <- result$correlation_abs
    result$.ymin <- result$conf.low_abs
    result$.ymax <- result$conf.high_abs
    y_lab <- sprintf("Absolute %s correlation with %s",
                     tools::toTitleCase(method), formatter(response))
    y_breaks <- seq(0, 1, by = 0.1)
  } else {
    result$.y    <- result$correlation
    result$.ymin <- result$conf.low
    result$.ymax <- result$conf.high
    y_lab <- sprintf("%s correlation with %s",
                     tools::toTitleCase(method), formatter(response))
    y_breaks <- ggplot2::waiver()
  }

  # --- Display labels (R4) ------------------------------------------------
  # Metric stem for x axis: format_metric_name strips the region suffix.
  # mrna-suffixed whole-transcript columns lose their " mRNA" tail here
  # (e.g. "Mean intron length", "uORF count"); standalone scalars have no
  # region token in the variable name so they pass through unchanged.
  result$metric_display <- format_metric_name(result$variable)

  # Within-facet ordering: by max |correlation| across regions, per supergroup,
  # per species. Apply globally then let facet_grid(scales="free_x") slice it.
  # Identifier guarantees uniqueness when the same metric_display appears in
  # different supergroups (rare but possible).
  result <- result |>
    dplyr::mutate(metric_key = paste(supergroup, metric_display, sep = "::"))

  order_tbl <- result |>
    dplyr::group_by(metric_key) |>
    dplyr::summarise(max_r = max(correlation_abs, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_r))

  result <- result |>
    dplyr::mutate(
      metric_key = factor(metric_key, levels = order_tbl$metric_key),
      region_f   = factor(region,
                          levels = intersect(names(region_colours),
                                             unique(region)))
    )

  # --- Plot ---------------------------------------------------------------
  dodge_width <- 0.5
  dodge       <- ggplot2::position_dodge(width = dodge_width)

  # Axis tick labels: strip the supergroup prefix that we added to keep the
  # factor levels unique. Lookup table for the scale.
  x_label_lookup <- setNames(
    sub("^[^:]+::", "", levels(result$metric_key)),
    levels(result$metric_key)
  )

  p <- ggplot2::ggplot(
    result,
    ggplot2::aes(
      x = metric_key, y = .y,
      colour = region_f, shape = region_f, fill = region_f
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .ymin, ymax = .ymax),
      width = 0.2, position = dodge, linewidth = 0.4
    ) +
    ggplot2::geom_point(size = 2.5, position = dodge, stroke = 0.5) +
    ggplot2::scale_colour_manual(
      values = region_colours, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r],
      name   = "Region"
    ) +
    ggplot2::scale_fill_manual(
      values = region_colours, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r],
      name   = "Region"
    ) +
    ggplot2::scale_shape_manual(
      values = region_shapes, drop = TRUE,
      labels = function(r) REGION_DISPLAYS[r],
      name   = "Region"
    ) +
    ggplot2::scale_x_discrete(labels = x_label_lookup) +
    ggplot2::scale_y_continuous(breaks = y_breaks) +
    ggplot2::labs(
      title    = sprintf("%s correlation with %s",
                         tools::toTitleCase(method), formatter(response)),
      subtitle = sprintf(
        paste0("%s with %d%% CI; ",
               "within-facet order by max |r|"),
        tools::toTitleCase(method), round(conf * 100)
      ),
      x = NULL,
      y = y_lab
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title        = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle     = ggplot2::element_text(size = 10, colour = "grey30"),
      axis.text.x       = ggplot2::element_text(angle = 315, vjust = 0.5,
                                                hjust = 0, size = 8),
      panel.background  = ggplot2::element_rect(fill = "grey95"),
      panel.grid.major.y = ggplot2::element_line(colour = "grey75",
                                                 linetype = "dotted"),
      panel.grid.minor.y = ggplot2::element_line(colour = "grey90",
                                                 linetype = "dotted"),
      panel.grid.major.x = ggplot2::element_line(colour = "grey90",
                                                 linetype = "dashed"),
      strip.background  = ggplot2::element_rect(fill = "grey90",
                                                colour = "black"),
      strip.text        = ggplot2::element_text(face = "bold")
    )

  # Reference line and zero line
  if (!is.null(sig_threshold)) {
    p <- p + ggplot2::geom_hline(
      yintercept = sig_threshold,
      linetype = "dotted", colour = "red", linewidth = 0.4
    )
  }
  if (!absolute) {
    p <- p + ggplot2::geom_hline(
      yintercept = 0, linetype = "solid", colour = "grey40", linewidth = 0.3
    )
  }

  # --- Faceting -----------------------------------------------------------
  if (has_species) {
    p <- p + ggplot2::facet_grid(
      species ~ supergroup,
      scales = "free_x", space = "free_x"
    )
  } else {
    p <- p + ggplot2::facet_grid(
      ~ supergroup,
      scales = "free_x", space = "free_x"
    )
  }

  # --- R9: return clean table ---------------------------------------------
  table_out <- result |>
    dplyr::select(dplyr::any_of(c(
      "species", "variable", "group", "supergroup",
      "metric_stem", "metric_display", "region",
      "n", "correlation", "conf.low", "conf.high",
      "p_value", "q_value",
      "correlation_abs", "conf.low_abs", "conf.high_abs"
    )))

  list(plot = p, table = table_out)
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

  # Reusable runner: response, suffix for filenames, sig threshold
  # Each job: response variable, filename suffix, sig threshold, group filter,
  # top-N filter, output width (mm). The codon/AA job overrides groups to opt
  # IN to the high-cardinality groups that the default expansion excludes, and
  # uses top_n_per_group to keep the figure legible.
  # Reusable narrowings ("the two reported NMD columns", "the four core
  # length columns") live as bundles in GROUP_BUNDLES (config.R), not inline
  # here. The broad jobs pull the full `intrinsic`/`splicing` supergroups and
  # additionally name the bundles so the bundles' pick lists attach to the
  # `lengths`/`nmd` group keys those supergroups bring in.
  broad_groups <- c("structure", "intrinsic", "splicing", "regulatory",
                    "codon_freqs", "aa_freqs", "nuc_ratios",
                    "nmd_reported", "lengths_core")

  jobs <- list(
    list(response = "halflife",
         suffix   = "halflife",
         sig      = 0.02,
         groups   = broad_groups,
         absolute = TRUE,
         top_n    = list(codon_freqs = 3, aa_freqs = 3),
         width    = 380),
    list(response = "translation_efficiency",
         suffix   = "translation_efficiency",
         sig      = 0.02,
         groups   = broad_groups,
         absolute = TRUE,
         top_n    = list(codon_freqs = 3, aa_freqs = 3),
         width    = 380),
    list(response = "halflife",
         suffix   = "halflife_codon_aa",
         sig      = 0.02,
         groups   = c("codon_freqs", "aa_freqs"),
         absolute = TRUE,
         top_n    = list(codon_freqs = 15, aa_freqs = 10),
         width    = 200),
    list(response = "halflife",
         suffix   = "halflife_nuc_ratios",
         sig      = 0.02,
         groups   = c("nuc_ratios"),
         absolute = TRUE,
         top_n    = list(),
         width    = 200)
  )

  for (job in jobs) {
    if (!job$response %in% names(df)) {
      message("Skipping: ", job$response, " not in dataset")
      next
    }

    message("\nDotplot for: ", job$response, " (", job$suffix, ")")
    out <- feature_correlation_dotplot(
      df,
      response        = job$response,
      sig_threshold   = job$sig,
      groups          = job$groups,
      absolute        = job$absolute,
      top_n_per_group = job$top_n
    )
    print(out$plot)
    
    ggplot2::ggsave(
      file.path(OUTPUT_DIR, "plots",
                paste0("feature_correlation_dotplot_", job$suffix, ".jpg")),
      plot = out$plot,
      width = job$width, height = 200, units = "mm", dpi = 300
    )
    write.csv(
      out$table,
      file.path(OUTPUT_DIR, "tables",
                paste0("feature_correlation_dotplot_", job$suffix, ".csv")),
      row.names = FALSE
    )
  }

  message("\nDotplots complete:")
  message("  ", file.path(OUTPUT_DIR, "plots"),
          "/feature_correlation_dotplot_*.jpg")
  message("  ", file.path(OUTPUT_DIR, "tables"),
          "/feature_correlation_dotplot_*.csv")
}
