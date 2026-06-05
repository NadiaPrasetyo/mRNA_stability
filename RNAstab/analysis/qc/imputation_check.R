# =============================================================================
# QC: mRNA MFE imputation diagnostics
# =============================================================================
# Checks how well the length-weighted region-based imputation of mRNA z-score
# (and additive region-based imputation of mRNA score) agrees with the
# directly-computed mRNA-level values, on transcripts that have both.
#
# Fixes vs. old pipeline:
#   - The old impute_Zscores.R filtered out rows with any NA region z-score
#     before imputing, which defeats the purpose of imputation. The corrected
#     pipeline imputes from whatever regions are present.
#   - Labels here say "Spearman \u03c1" not "Pearson" — the old code used
#     Spearman but labelled plots as Pearson.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

species <- "human"
df <- build_dataset(species)

# The imputation step writes both the raw mRNA z-score (from direct folding
# of the full mRNA) and an *_imputed variant (the length-weighted region avg).
# The coalesce step has already filled rnafold_zscore_mrna with the imputed value
# where the direct one was missing, but the imputed column is preserved so
# we can compare them here on transcripts where BOTH exist.

need <- c("rnafold_zscore_mrna", "rnafold_zscore_mrna_imputed",
          "rnafold_score_mrna",  "rnafold_score_mrna_imputed")
missing <- setdiff(need, names(df))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "),
       ". Rebuild the dataset with build_dataset('", species, "', rebuild = TRUE).")
}

cmp <- df |>
  filter(!is.na(rnafold_zscore_mrna), !is.na(rnafold_zscore_mrna_imputed),
         !is.na(rnafold_score_mrna),  !is.na(rnafold_score_mrna_imputed))

rho_z  <- cor(cmp$rnafold_zscore_mrna, cmp$rnafold_zscore_mrna_imputed, method = "spearman")
rho_sc <- cor(cmp$rnafold_score_mrna,  cmp$rnafold_score_mrna_imputed,  method = "spearman")


# --- Plot 1: z-score — direct vs imputed -------------------------------------
p_zscore <- ggplot(cmp, aes(x = rnafold_zscore_mrna, y = rnafold_zscore_mrna_imputed)) +
  geom_point(alpha = 0.7, colour = "steelblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "firebrick", linewidth = 1) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Spearman \u03c1 = %.3f", rho_z),
           hjust = -0.1, vjust = 1.5, size = 5) +
  labs(
    title = "mRNA MFE z-score — direct vs length-weighted imputation",
    x = "Direct (mRNA folding)",
    y = "Imputed (length-weighted region average)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )


# --- Plot 2: score — direct vs imputed ---------------------------------------
p_score <- ggplot(cmp, aes(x = rnafold_score_mrna, y = rnafold_score_mrna_imputed)) +
  geom_point(alpha = 0.7, colour = "steelblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "firebrick", linewidth = 1) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Spearman \u03c1 = %.3f", rho_sc),
           hjust = -0.1, vjust = 1.5, size = 5) +
  labs(
    title = "mRNA MFE score — direct vs additive imputation",
    x = "Direct (mRNA folding)",
    y = "Imputed (sum of region scores)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )


# --- Plot 3: region MFE delta distributions ---------------------------------
# Distribution of (observed - expected) by region, per the thermodynamic model.
delta_cols <- grep("^mfe_delta_", names(df), value = TRUE)
delta_cols <- setdiff(delta_cols, "mfe_delta_mrna")       # cleaner contrast

if (length(delta_cols) >= 1 && "mfe_delta_mrna" %in% names(df)) {
  long_delta <- df |>
    select(all_of(c(delta_cols, "mfe_delta_mrna"))) |>
    pivot_longer(everything(),
                 names_to  = "region",
                 values_to = "delta") |>
    mutate(region = sub("^mfe_delta_", "", region)) |>
    filter(!is.na(delta))

  p_delta <- ggplot(long_delta, aes(x = delta, fill = region)) +
    geom_density(alpha = 0.5, colour = "black") +
    coord_cartesian(xlim = c(-60, 35)) +
    labs(
      title = "MFE delta (observed - expected) by region",
      x = "MFE delta",
      y = "Density",
      fill = "Region"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text  = element_text(size = 12)
    )
  print(p_delta)
}

print(p_zscore)
print(p_score)
