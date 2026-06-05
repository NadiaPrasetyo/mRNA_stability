# =============================================================================
# QC: Observed vs expected MFE
# =============================================================================
# Diagnostic plots for the Agarwal & Kelley thermodynamic model:
#   Expected MFE = (a · GC^b + c) · length + d
#
# If the model is well-calibrated, observed and expected should lie near the
# y = x line. These plots also colour by GC content and length to expose any
# remaining residual structure.
#
# Run this *after* build_dataset() has populated the cache.
# =============================================================================

source("R/load_all.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

species <- "human"                              # change as needed
region  <- "mrna"                               # 5utr / cds / 3utr / mrna / con53
df <- build_dataset(species)

score_col <- paste0("rnafold_score_",    region)
exp_col   <- paste0("mfe_expected_", region)
gc_col    <- paste0("gc_",           region)
len_col   <- paste0("length_",       region)
delta_col <- paste0("mfe_delta_",    region)

req <- c(score_col, exp_col, gc_col, len_col)
stopifnot(all(req %in% names(df)))

complete <- df |>
  filter(!is.na(.data[[score_col]]), !is.na(.data[[exp_col]]))

rho <- cor(complete[[score_col]], complete[[exp_col]], method = "spearman")
rho_label <- sprintf("Spearman \u03c1 = %.3f", rho)

x_lims <- quantile(complete[[score_col]], probs = c(0.01, 0.99), na.rm = TRUE)
y_lims <- quantile(complete[[exp_col]],   probs = c(0.01, 0.99), na.rm = TRUE)

# --- Plot 1: coloured by GC ---------------------------------------------------
plot_by_gc <- ggplot(complete,
                     aes(x = .data[[score_col]], y = .data[[exp_col]])) +
  geom_point(aes(colour = .data[[gc_col]]), alpha = 0.8, size = 1) +
  scale_colour_viridis_c(option = "plasma") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "red", linewidth = 1) +
  annotate("text", x = -Inf, y = Inf, label = rho_label,
           hjust = -0.2, vjust = 2, size = 4) +
  coord_cartesian(xlim = x_lims, ylim = y_lims) +
  labs(
    title    = paste0("Observed vs modelled MFE — ", format_col_name(region)),
    subtitle = bquote(Modelled ~ MFE == (-0.84 ~ GC^2.35 - 0.15) %.% length + 13.56),
    x = format_col_name(score_col),
    y = format_col_name(exp_col),
    colour = format_col_name(gc_col)
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title    = element_text(size = 16),
    axis.text     = element_text(size = 14)
  )

# --- Plot 2: coloured by length -----------------------------------------------
plot_by_length <- plot_by_gc %+%
  complete +
  aes(colour = .data[[len_col]])
plot_by_length <- plot_by_length +
  scale_colour_viridis_c(option = "viridis") +
  labs(colour = format_col_name(len_col))

# --- Plot 3: delta vs length, coloured by GC ---------------------------------
if (delta_col %in% names(complete)) {
  plot_delta_vs_length <- ggplot(complete,
                                 aes(x = .data[[len_col]],
                                     y = .data[[delta_col]])) +
    geom_point(aes(colour = .data[[gc_col]]), alpha = 0.7, size = 2) +
    scale_colour_viridis_c(option = "viridis") +
    labs(
      title  = paste0("Length vs MFE \u0394 — ", format_col_name(region)),
      x      = format_col_name(len_col),
      y      = format_col_name(delta_col),
      colour = format_col_name(gc_col)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14)
    )
}

print(plot_by_gc)
print(plot_by_length)
if (exists("plot_delta_vs_length")) print(plot_delta_vs_length)
