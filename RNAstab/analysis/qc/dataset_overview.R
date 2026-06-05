# =============================================================================
# Dataset QC overview
# =============================================================================
# Four diagnostic views of the built dataset, run before any downstream
# analysis. Each function is independently usable and returns list(plot, table)
# per R9 of PIPELINE_GUIDE.md.
#
#   1. halflife_distribution_plot   — shape of the response, per species
#   2. missingness_by_group_plot    — coverage heatmap across FEATURE_PATTERNS
#   3. sample_size_summary_plot     — n transcripts at each completeness stage
#   4. loader_health_plot           — canary column coverage per loader
#
# All functions facet/dodge on `species` (R6), so they are already correct for
# multi-species inputs from build_all().
#
# Usage:
#   source("R/load_all.R")
#   source("analysis/qc/dataset_overview.R")
#   df  <- build_dataset("human")
#   out <- halflife_distribution_plot(df); print(out$plot)
#
# Or run end-to-end:
#   Rscript analysis/qc/dataset_overview.R
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(purrr)
  library(tibble)
})


# -----------------------------------------------------------------------------
# 1. Halflife distribution
# -----------------------------------------------------------------------------

#' Distribution of halflife, faceted by species.
#'
#' Halflife from the Agarwal paper is PC1 of quantile-normalised log half-lives,
#' so the natural QC question is whether it is approximately symmetric and
#' centred — not whether to log-transform.
#'
#' @param df         Dataframe from build_dataset() or build_all().
#' @param formatter  Display formatter (default format_col_name).
#' @return list(plot, table). Table: species, n, mean, median, sd, q25, q75,
#'   min, max, skew.
halflife_distribution_plot <- function(df, formatter = format_col_name) {

  # R5: guard
  if (!"halflife" %in% names(df)) {
    stop("halflife missing — wrong species, or response variable not joined")
  }
  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }

  d <- df |> dplyr::filter(!is.na(halflife))
  if (nrow(d) == 0) stop("No non-NA halflife rows")

  summary_tbl <- d |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      n      = dplyr::n(),
      mean   = mean(halflife),
      median = stats::median(halflife),
      sd     = stats::sd(halflife),
      q25    = stats::quantile(halflife, 0.25),
      q75    = stats::quantile(halflife, 0.75),
      min    = min(halflife),
      max    = max(halflife),
      skew   = mean((halflife - mean(halflife))^3) / stats::sd(halflife)^3,
      .groups = "drop"
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(x = halflife)) +
    ggplot2::geom_density(fill = "#4B0082", alpha = 0.3, colour = "#4B0082") +
    ggplot2::geom_rug(alpha = 0.05, length = ggplot2::unit(0.015, "npc")) +
    ggplot2::geom_vline(
      data = summary_tbl,
      ggplot2::aes(xintercept = median),
      linetype = "dashed", colour = "#4B0082"
    ) +
    ggplot2::geom_text(
      data = summary_tbl,
      ggplot2::aes(
        x = -Inf, y = Inf,
        label = sprintf("n = %s\nmedian = %.2f\nskew = %.2f",
                        scales::comma(n), median, skew)
      ),
      hjust = -0.1, vjust = 1.2, size = 3, inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~ species, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title    = paste0("Distribution of ", formatter("halflife")),
      subtitle = "Agarwal aggregated: PC1 of quantile-normalised log half-lives — expected approx. symmetric, centred near zero",
      x        = formatter("halflife"),
      y        = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 10),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  list(plot = p, table = summary_tbl)
}


# -----------------------------------------------------------------------------
# 2. Missingness by feature group
# -----------------------------------------------------------------------------

#' Tile heatmap of mean % non-NA per feature group × species.
#'
#' For each group in FEATURE_PATTERNS, computes the average non-NA rate
#' across all member columns, per species. Tile text shows percentage and
#' the number of columns in the group.
#'
#' @param df         Dataframe from build_dataset() or build_all().
#' @param groups     Character vector of FEATURE_PATTERNS keys. Defaults to
#'                   all groups except `lengths_some` (redundant subset of
#'                   `lengths`).
#' @param formatter  Display formatter (default format_col_name).
#' @return list(plot, table). Table: species, group, n_columns,
#'   mean_nonna_pct, min_nonna_pct, n_columns_all_na.
missingness_by_group_plot <- function(df,
                                      groups    = NULL,
                                      formatter = format_col_name) {

  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }

  if (is.null(groups)) {
    # Exclude lengths_some (redundant subset of lengths). rnafold_scores and
    # mfe_scores are aliases; same for rnafold_zscores / mfe_zscores.
    groups <- setdiff(names(FEATURE_PATTERNS),
                      c("lengths_some", "mfe_scores", "mfe_zscores"))
  }

  result <- purrr::map_dfr(unique(df$species), function(sp) {
    sub <- df |> dplyr::filter(species == sp)
    purrr::map_dfr(groups, function(g) {

      # R3: use fg_columns, not hand-rolled regex
      cols <- fg_columns(sub, g)

      if (length(cols) == 0) {
        tibble::tibble(
          species         = sp,
          group           = g,
          n_columns       = 0L,
          mean_nonna_pct  = NA_real_,
          min_nonna_pct   = NA_real_,
          n_columns_all_na = 0L
        )
      } else {
        nonna_pct <- vapply(cols, function(co) {
          100 * mean(!is.na(sub[[co]]))
        }, numeric(1))

        tibble::tibble(
          species          = sp,
          group            = g,
          n_columns        = length(cols),
          mean_nonna_pct   = mean(nonna_pct),
          min_nonna_pct    = min(nonna_pct),
          n_columns_all_na = sum(nonna_pct == 0)
        )
      }
    })
  })

  # Order y-axis by overall coverage (worst at bottom)
  result <- result |>
    dplyr::mutate(
      # 1. Convert to factor
      group = factor(group),
      # 2. Drop levels that have no columns (where n_columns == 0)
      # Or just drop levels that aren't present in this specific result set
      group = forcats::fct_drop(group),
      # 3. Now reorder safely
      group = forcats::fct_reorder(group, mean_nonna_pct, 
                                   .fun = mean, .na_rm = TRUE)
    )
  p <- ggplot2::ggplot(
    result,
    ggplot2::aes(x = species, y = group, fill = mean_nonna_pct)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(is.na(mean_nonna_pct), "no cols",
                                  sprintf("%.0f%%\n(%d col%s)",
                                          mean_nonna_pct, n_columns,
                                          ifelse(n_columns == 1, "", "s")))),
      size = 3, colour = "black"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#b2182b", mid = "#f7f7f7", high = "#2166ac",
      midpoint = 50, limits = c(0, 100),
      name = "% non-NA",
      na.value = "grey90"
    ) +
    ggplot2::labs(
      title    = "Feature-group coverage",
      subtitle = "Mean % non-NA across columns in each group; cells with 0 cols mean no FEATURE_PATTERNS match in this build",
      x        = formatter("species"),
      y        = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 10),
      axis.text.y      = ggplot2::element_text(size = 10),
      panel.grid       = ggplot2::element_blank()
    )

  list(plot = p, table = result)
}


# -----------------------------------------------------------------------------
# 3. Sample size summary
# -----------------------------------------------------------------------------

#' Dodged bar chart of transcript counts at each completeness stage.
#'
#' Stages are cumulative-ish filters: total, with halflife, with CAI, with TE,
#' with mRNA MFE z-score, and the intersection of all of those (a Saluki-style
#' core feature set).
#'
#' @param df         Dataframe from build_dataset() or build_all().
#' @param formatter  Display formatter (default format_col_name).
#' @return list(plot, table). Table: one row per species, one column per stage.
sample_size_summary_plot <- function(df, formatter = format_col_name) {

  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }

  # R5: guard each column individually
  has_hl  <- "halflife"               %in% names(df)
  has_cai <- "cai"                    %in% names(df)
  has_te  <- "translation_efficiency" %in% names(df)
  has_mfe <- "rnafold_zscore_mrna"    %in% names(df)

  summary_tbl <- df |>
    dplyr::group_by(species) |>
    dplyr::summarise(
      total           = dplyr::n(),
      with_halflife   = if (has_hl)  sum(!is.na(halflife))                else NA_integer_,
      with_cai        = if (has_cai) sum(!is.na(cai))                     else NA_integer_,
      with_te         = if (has_te)  sum(!is.na(translation_efficiency))  else NA_integer_,
      with_mrna_mfe_z = if (has_mfe) sum(!is.na(rnafold_zscore_mrna))     else NA_integer_,
      complete_core   = if (has_hl && has_cai && has_te && has_mfe) {
        sum(!is.na(halflife) & !is.na(cai) &
            !is.na(translation_efficiency) & !is.na(rnafold_zscore_mrna))
      } else NA_integer_,
      .groups = "drop"
    )

  stage_levels <- c("total", "with_halflife", "with_cai",
                    "with_te", "with_mrna_mfe_z", "complete_core")
  stage_labels <- c(
    total           = "All transcripts",
    with_halflife   = "With half-life",
    with_cai        = "With CAI",
    with_te         = "With TE",
    with_mrna_mfe_z = "With mRNA MFE z",
    complete_core   = "Complete (HL+CAI+TE+MFE)"
  )

  long <- summary_tbl |>
    tidyr::pivot_longer(-species, names_to = "stage", values_to = "n") |>
    dplyr::filter(!is.na(n)) |>
    dplyr::mutate(
      stage       = factor(stage, levels = stage_levels),
      stage_label = factor(stage_labels[as.character(stage)],
                           levels = unname(stage_labels))
    )

  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = stage_label, y = n, fill = species)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::geom_text(
      ggplot2::aes(label = scales::comma(n)),
      position = ggplot2::position_dodge(width = 0.9),
      vjust = -0.3, size = 3
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = formatter("species")) +
    ggplot2::scale_y_continuous(labels = scales::comma_format(),
                                expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(
      title    = "Sample size by completeness stage",
      subtitle = "Each bar is the count of transcripts with non-NA value(s) for the named filter",
      x        = NULL,
      y        = "Transcripts"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 10),
      axis.text.x        = ggplot2::element_text(angle = 25, hjust = 1),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )

  list(plot = p, table = summary_tbl)
}


# -----------------------------------------------------------------------------
# 4. Loader health (canary columns)
# -----------------------------------------------------------------------------

#' Bar chart of % non-NA for one representative ("canary") column per loader.
#'
#' Loaders skip missing input files silently, returning NULL — so an entire
#' loader can silently fail to contribute. This plot picks one column per
#' loader and shows its coverage. 0% means the loader's input file is missing
#' or the join key failed.
#'
#' To add or change canaries: edit the `canaries` vector below.
#'
#' @param df         Dataframe from build_dataset() or build_all().
#' @param formatter  Display formatter (default format_col_name).
#' @return list(plot, table). Table: species, loader, canary, present, nonna_pct.
loader_health_plot <- function(df, formatter = format_col_name) {

  if (!"species" %in% names(df)) {
    stop("species column missing — pipeline invariant violated")
  }

  # Map: human-readable loader label -> a column it should produce.
  # Picked one column per loader from the schema in PIPELINE_GUIDE §2.4.
  canaries <- c(
    "RNAfold (global MFE)"        = "rnafold_score_mrna",
    "RNAfold z-scores"            = "rnafold_zscore_cds",
    "RNALfold (local MFE)"        = "rnalfold_score_cds",
    "RNAup (UTR pair)"            = "rnaup_score_utrpair",
    "Sequence basic (GC/length)"  = "gc_content_cds",
    "Stop-free regions"           = "stopfree_cds",
    "NMD"                         = "nmd_snv_fragile_codon_density_mrna",
    "Architecture (introns)"      = "intron_mean",
    "Junctions / distances"       = "n_cds_junctions",
    "uORFs"                       = "n_uorfs",
    "Amino-acid frequencies"      = "aa_freq_leu",
    "CAI"                         = "cai",
    "Translation efficiency"      = "translation_efficiency",
    "Agarwal features"            = "expression",
    "Saluki predictions"          = "saluki_prediction",
    "Dani probing (icSHAPE)"      = "shape_mean_5utr",
    "Dani probing (Keth-seq)"     = "keth_mean_5utr"
  )

  result <- purrr::map_dfr(unique(df$species), function(sp) {
    sub <- df |> dplyr::filter(species == sp)
    purrr::imap_dfr(canaries, function(col, lab) {
      # R5: guard column existence before access
      if (!col %in% names(sub)) {
        tibble::tibble(
          species   = sp,
          loader    = lab,
          canary    = col,
          present   = FALSE,
          nonna_pct = 0
        )
      } else {
        tibble::tibble(
          species   = sp,
          loader    = lab,
          canary    = col,
          present   = TRUE,
          nonna_pct = 100 * mean(!is.na(sub[[col]]))
        )
      }
    })
  })

  # Order loaders by max coverage across species (best at top)
  result <- result |>
    dplyr::mutate(loader = forcats::fct_reorder(loader, nonna_pct, .fun = max))

  p <- ggplot2::ggplot(
    result,
    ggplot2::aes(x = nonna_pct, y = loader, fill = species)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.0f%%", nonna_pct)),
      position = ggplot2::position_dodge(width = 0.7),
      hjust = -0.15, size = 3
    ) +
    ggplot2::scale_fill_brewer(palette = "Set2", name = formatter("species")) +
    ggplot2::scale_x_continuous(
      limits = c(0, 115),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%")
    ) +
    ggplot2::labs(
      title    = "Loader health: canary column coverage",
      subtitle = "% of transcripts with non-NA value for one representative column per loader; 0% = file missing or join key failed",
      x        = "Transcripts with non-NA canary value",
      y        = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title         = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle      = ggplot2::element_text(size = 10),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )

  list(plot = p, table = result)
}


# -----------------------------------------------------------------------------
# Runner: top-to-bottom execution writes all four plots + tables
# -----------------------------------------------------------------------------

if (sys.nframe() == 0 || identical(environment(), globalenv())) {

  # Switch to build_all() once mouse (or others) come online.
  df <- build_dataset("human")

  # R8: outputs go under OUTPUT_DIR
  dir.create(file.path(OUTPUT_DIR, "plots"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "tables"),
             showWarnings = FALSE, recursive = TRUE)

  # Each entry: function + dimensions tuned to its content
  jobs <- list(
    halflife_distribution = list(
      fn = halflife_distribution_plot, w = 210, h = 148
    ),
    missingness_by_group  = list(
      fn = missingness_by_group_plot,  w = 210, h = 250
    ),
    sample_size_summary   = list(
      fn = sample_size_summary_plot,   w = 210, h = 148
    ),
    loader_health         = list(
      fn = loader_health_plot,         w = 210, h = 230
    )
  )

  for (nm in names(jobs)) {
    message("QC: ", nm)
    out <- jobs[[nm]]$fn(df)
    print(out$plot)

    ggplot2::ggsave(
      file.path(OUTPUT_DIR, "plots", paste0("qc_", nm, ".jpg")),
      plot = out$plot,
      width = jobs[[nm]]$w, height = jobs[[nm]]$h,
      units = "mm", dpi = 300
    )
    write.csv(
      out$table,
      file.path(OUTPUT_DIR, "tables", paste0("qc_", nm, ".csv")),
      row.names = FALSE
    )
  }

  message("\nQC complete. Outputs:")
  message("  ", file.path(OUTPUT_DIR, "plots"),  "/qc_*.jpg")
  message("  ", file.path(OUTPUT_DIR, "tables"), "/qc_*.csv")
}
