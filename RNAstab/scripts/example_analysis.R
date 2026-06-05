# =============================================================================
# Example analysis
# =============================================================================
# A minimal worked example showing how to consume the pipeline:
#   1. Load datasets (cached → instant on re-runs).
#   2. Stack species for cross-species comparisons.
#   3. Use feature-group helpers to grab the right columns.
#   4. Plot / model.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})


# --- 1. Single species ------------------------------------------------------
human <- build_dataset("human")
cat("Human:", nrow(human), "transcripts,", ncol(human), "features\n")


# --- 2. Cross-species -------------------------------------------------------
# bind_rows stacks them on a shared `species` column. Columns absent from
# one species are NA for those rows.
combined <- build_all()

# Compare half-life distributions across species
combined |>
  filter(!is.na(halflife)) |>
  group_by(species) |>
  summarise(
    n         = n(),
    median_hl = median(halflife),
    iqr_hl    = IQR(halflife)
  ) |>
  print()


# --- 3. Feature groups in action --------------------------------------------
# Pull every MFE z-score for human:
mfe_z <- human |> select(transcript_id, fg("rnafold_zscores"))
cat("MFE z-score columns:\n"); print(names(mfe_z))

# Combine several groups into a modelling frame:
model_df <- human |>
  select(
    halflife,
    fg("rnafold_zscores"),
    fg("mfe_deltas"),
    fg("rnalfold_zscores")
  ) |>
  na.omit()

cat("\nModel frame:", nrow(model_df), "rows,", ncol(model_df), "columns\n")


# --- 4. Quick scatter -------------------------------------------------------
if ("rnafold_zscore_mrna" %in% names(human) && "halflife" %in% names(human)) {
  p <- ggplot(human, aes(x = rnafold_zscore_mrna, y = halflife)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, colour = "firebrick") +
    labs(
      x     = format_col_name("rnafold_zscore_mrna"),
      y     = format_col_name("halflife"),
      title = "mRNA MFE z-score vs half-life (human)"
    ) +
    theme_minimal()

  ggsave(file.path(OUTPUT_DIR, "plots", "example_mfe_vs_halflife.jpg"),
         plot = p, width = 210, height = 148, units = "mm", dpi = 300)
  cat("Saved example plot to", file.path(OUTPUT_DIR, "plots"), "\n")
}

