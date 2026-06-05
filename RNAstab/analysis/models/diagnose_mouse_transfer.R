# diagnose_mouse_transfer.R
#
# Cross-species transfer diagnostics: apply the human-trained half-life model
# to a mouse dataset and produce parallel plots to diagnose_halflife_model.R.
# No refitting and no split — the saved final_full_fit (trained on all human
# data) is used directly to predict on mouse.
#
#   Dashboard (2x3, saved as mouse_transfer_dashboard.{png,pdf}):
#     1. Predicted vs observed (hexbin)
#     2. Residuals vs predicted
#     3. Residual distribution (mouse, with human held-out as reference)
#     4. Calibration by predicted decile
#     5. Observed half-life within each predicted decile (violin + box)
#     6. Top-K recovery curve
#
#   Combined filtering view:
#     filtering_quality_mouse.{png,pdf}: top-K + absolute-cutoff side-by-side
#
# Also saves mouse predictions to mouse_predictions.csv.
#
# Usage: edit the paths below, then: Rscript diagnose_mouse_transfer.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(bonsai)
  library(patchwork)
})


# ----------------------------- Configuration --------------------------------

MODEL_PATH <- "models/halflife_lgbm_v5.rds"
MOUSE_PATH <- "data/cache/mouse_dataset_v5.rds"
OUT_DIR    <- "figures/mouse_transfer"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ----------------------------- 1. Load and schema check ---------------------

m     <- readRDS(MODEL_PATH)
mouse <- readRDS(MOUSE_PATH)

target_col   <- m$target_col
feature_cols <- m$feature_cols

missing <- setdiff(feature_cols, names(mouse))
if (length(missing) > 0) {
  cat("\nMouse data is missing", length(missing),
      "features required by the model:\n")
  print(missing)
  stop("Cannot proceed - regenerate mouse features upstream first.")
}

if (!target_col %in% names(mouse)) {
  stop(sprintf("Mouse data is missing the target column '%s'.", target_col))
}

cat("Schema OK:\n")
cat("  features required: ", length(feature_cols), "\n")
cat("  rows in mouse:     ", nrow(mouse), "\n")


# ----------------------------- 2. Predict on mouse --------------------------
# m$workflow is the deployment fit (final_wf fit on all human data).

mouse_eval <- mouse %>% filter(!is.na(.data[[target_col]]))
cat("  rows with measured half-life: ", nrow(mouse_eval), "\n\n")

cat("Predicting on mouse with the human-trained model...\n")
mouse_pred <- augment(m$workflow, mouse_eval) %>%
  filter(!is.na(.pred), !is.na(.data[[target_col]])) %>%
  mutate(residual = .pred - .data[[target_col]])


# ----------------------------- 3. Summary statistics ------------------------

stats <- function(d, truth_col, pred_col = ".pred") {
  truth <- d[[truth_col]]; pred <- d[[pred_col]]
  list(
    rmse = sqrt(mean((truth - pred)^2)),
    mae  = mean(abs(truth - pred)),
    r    = cor(truth, pred),
    rsq  = cor(truth, pred)^2,
    rho  = cor(truth, pred, method = "spearman"),
    n    = length(truth)
  )
}

s_mouse <- stats(mouse_pred, target_col)

# Reference: human held-out (from saved model)
saved_rsq  <- m$test_metrics %>% filter(.metric == "rsq")  %>% pull(.estimate)
saved_rmse <- m$test_metrics %>% filter(.metric == "rmse") %>% pull(.estimate)
saved_rho  <- m$test_spearman

cat(sprintf("\nMouse transfer:           n=%d, R^2=%.3f, Spearman=%.3f, RMSE=%.3f\n",
            s_mouse$n, s_mouse$rsq, s_mouse$rho, s_mouse$rmse))
cat(sprintf("Human held-out reference:        R^2=%.3f, Spearman=%.3f, RMSE=%.3f\n",
            saved_rsq, saved_rho, saved_rmse))
cat(sprintf("Transfer change:                 R^2 %+0.3f, Spearman %+0.3f, RMSE %+0.3f\n\n",
            s_mouse$rsq - saved_rsq,
            s_mouse$rho - saved_rho,
            s_mouse$rmse - saved_rmse))


# ----------------------------- 4. Plot helpers ------------------------------

theme_diag <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(colour = "grey40", size = 10),
      legend.position  = "right"
    )
}

stat_label <- function(s) {
  sprintf("R\u00B2 = %.3f   \u03C1 = %.3f   RMSE = %.2f   n = %d",
          s$rsq, s$rho, s$rmse, s$n)
}


# ----------------------------- 5. Plot 1: Predicted vs observed -------------

p_truth <- ggplot(mouse_pred, aes(x = .data[[target_col]], y = .pred)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  geom_abline(slope = 1, intercept = 0, colour = "white",
              linetype = "dashed", linewidth = 0.6) +
  geom_smooth(method = "lm", colour = "firebrick", se = FALSE,
              linewidth = 0.6, formula = y ~ x) +
  labs(title    = "Predicted vs observed (mouse transfer)",
       subtitle = stat_label(s_mouse),
       x        = "Observed mouse half-life",
       y        = "Predicted (from human-trained model)") +
  theme_diag()


# ----------------------------- 6. Plot 2: Residuals vs predicted ------------

p_resid <- ggplot(mouse_pred, aes(x = .pred, y = residual)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  geom_hline(yintercept = 0, colour = "white",
             linetype = "dashed", linewidth = 0.6) +
  geom_smooth(method = "loess", colour = "firebrick", se = FALSE,
              linewidth = 0.6, formula = y ~ x) +
  labs(title    = "Residuals vs predicted",
       subtitle = "Loess curve should hover around zero (slope = bias by predicted range)",
       x        = "Predicted",
       y        = "Residual (predicted \u2212 observed)") +
  theme_diag()


# ----------------------------- 7. Plot 3: Residual distribution -------------
# Single mouse density with a reference normal centred at 0 with sd = human
# held-out RMSE - the dashed curve is "what we'd see if mouse transferred
# as cleanly as the human held-out did". Comparing visually shows whether
# the mouse residuals are wider (transfer noisier) or shifted (systemic bias).

p_rdist <- ggplot(mouse_pred, aes(x = residual)) +
  geom_density(fill = "steelblue", alpha = 0.4, colour = "steelblue",
               linewidth = 0.6) +
  stat_function(
    fun    = dnorm,
    args   = list(mean = 0, sd = saved_rmse),
    colour = "grey40", linetype = "dashed", linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dotted",
             linewidth = 0.4) +
  geom_vline(xintercept = mean(mouse_pred$residual),
             colour = "firebrick", linewidth = 0.4) +
  labs(
    title    = "Residual distribution (mouse)",
    subtitle = sprintf(
      "Mouse: mean = %.3f (red), sd = %.3f  |  Dashed = normal(0, %.3f) human held-out reference",
      mean(mouse_pred$residual), sd(mouse_pred$residual), saved_rmse
    ),
    x = "Residual (predicted \u2212 observed)",
    y = "Density"
  ) +
  theme_diag()


# ----------------------------- 8. Plot 4: Calibration -----------------------

cal_data <- mouse_pred %>%
  mutate(bin = ntile(.pred, 10)) %>%
  group_by(bin) %>%
  summarise(
    pred_mean = mean(.pred),
    obs_mean  = mean(.data[[target_col]]),
    obs_se    = sd(.data[[target_col]]) / sqrt(n()),
    n         = n(),
    .groups   = "drop"
  )

p_cal <- ggplot(cal_data, aes(x = pred_mean, y = obs_mean)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60",
              linetype = "dashed") +
  geom_errorbar(aes(ymin = obs_mean - obs_se,
                    ymax = obs_mean + obs_se),
                width = 0, colour = "grey50") +
  geom_point(aes(size = n), colour = "steelblue") +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(title    = "Calibration by predicted decile",
       subtitle = "Points should lie on the dashed line",
       x        = "Mean predicted in bin",
       y        = "Mean observed in bin") +
  theme_diag()


# ----------------------------- 9. Plot 5: Decile separability --------------

p_perfbin <- mouse_pred %>%
  mutate(bin = ntile(.pred, 10)) %>%
  ggplot(aes(x = factor(bin), y = .data[[target_col]])) +
  geom_violin(fill = "steelblue", alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(title    = "Observed mouse half-life within predicted decile",
       subtitle = "Less overlap between adjacent bins = better ranking",
       x        = "Predicted decile (1 = lowest, 10 = highest)",
       y        = "Observed mouse half-life") +
  theme_diag()


# ----------------------------- 10. Plot 6: Top-K recovery ------------------

ks <- seq(0.01, 1, by = 0.01)

recovery <- map_dfr(ks, function(k) {
  cutoff_pred  <- quantile(mouse_pred$.pred,         1 - k)
  cutoff_truth <- quantile(mouse_pred[[target_col]], 1 - k)
  predicted_top <- mouse_pred$.pred         >= cutoff_pred
  true_top      <- mouse_pred[[target_col]] >= cutoff_truth
  hits <- sum(predicted_top & true_top)
  tibble(
    k      = k,
    recall = hits / sum(true_top),
    chance = k
  )
})

p_topk <- ggplot(recovery, aes(x = k * 100)) +
  geom_line(aes(y = chance * 100, colour = "Random"),
            linetype = "dashed") +
  geom_line(aes(y = recall * 100, colour = "Model"), linewidth = 0.8) +
  scale_colour_manual(values = c(Model = "steelblue", Random = "grey50"),
                      name = NULL) +
  labs(title    = "Top-K recovery (filtering quality)",
       subtitle = "Top K% predicted: what % of the true top K% is captured?",
       x        = "Top K% by predicted half-life",
       y        = "% of true top K% recovered") +
  theme_diag() +
  coord_equal()


# ----------------------------- 11. Plot 7: Absolute-cutoff filtering ------

cutoffs <- seq(
  quantile(mouse_pred[[target_col]], 0.05, na.rm = TRUE),
  quantile(mouse_pred[[target_col]], 0.95, na.rm = TRUE),
  length.out = 60
)

threshold_data <- map_dfr(cutoffs, function(c) {
  predicted_pos <- mouse_pred$.pred         >= c
  true_pos      <- mouse_pred[[target_col]] >= c
  tp <- sum(predicted_pos & true_pos)
  tibble(
    cutoff     = c,
    precision  = if (sum(predicted_pos) > 0) tp / sum(predicted_pos) else NA_real_,
    recall     = if (sum(true_pos)      > 0) tp / sum(true_pos)      else NA_real_,
    prevalence = mean(true_pos)
  )
})

p_thresh <- threshold_data %>%
  pivot_longer(c(recall, precision, prevalence),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(
    metric,
    levels = c("recall", "precision", "prevalence"),
    labels = c("Recall (model)",
               "Precision (model)",
               "Prevalence (chance precision)")
  )) %>%
  ggplot(aes(x = cutoff, y = value * 100,
             colour = metric, linetype = metric)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(
    values = c("Recall (model)"                = "firebrick",
               "Precision (model)"             = "steelblue",
               "Prevalence (chance precision)" = "grey50"),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Recall (model)"                = "solid",
               "Precision (model)"             = "solid",
               "Prevalence (chance precision)" = "dashed"),
    name = NULL
  ) +
  labs(
    title    = "Filtering quality at absolute cutoffs (mouse)",
    subtitle = "Same cutoff applied to predicted and observed half-life score",
    x        = "Half-life score cutoff",
    y        = "%"
  ) +
  theme_diag() +
  theme(legend.position = "bottom")


# ----------------------------- 12. Combine into dashboard ------------------

dashboard <- (p_truth | p_resid | p_rdist) /
  (p_cal   | p_perfbin | p_topk) +
  plot_annotation(
    title    = "Half-life model - mouse transfer diagnostics",
    subtitle = sprintf(
      "Mouse n = %d  |  R\u00B2: mouse = %.3f vs human held-out = %.3f  (\u0394 = %+0.3f)",
      s_mouse$n, s_mouse$rsq, saved_rsq, s_mouse$rsq - saved_rsq
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(colour = "grey40")
    )
  )

ggsave(file.path(OUT_DIR, "mouse_transfer_dashboard.png"),
       dashboard, width = 16, height = 10, dpi = 200)
ggsave(file.path(OUT_DIR, "mouse_transfer_dashboard.pdf"),
       dashboard, width = 16, height = 10)


# ----------------------------- 13. Filtering quality combined --------------

filter_combined <- (p_topk | p_thresh) +
  plot_annotation(
    title    = "Filtering quality on mouse - relative vs absolute cutoffs",
    subtitle = "Same model, two framings: 'top K%' (rank-based) vs 'score \u2265 c' (threshold-based)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40")
    )
  )

ggsave(file.path(OUT_DIR, "filtering_quality_mouse.png"),
       filter_combined, width = 13, height = 5.5, dpi = 200)
ggsave(file.path(OUT_DIR, "filtering_quality_mouse.pdf"),
       filter_combined, width = 13, height = 5.5)


# ----------------------------- 14. Save predictions ------------------------

mouse_pred %>%
  select(any_of(m$meta_cols), all_of(target_col), .pred, residual) %>%
  write_csv(file.path(OUT_DIR, "mouse_predictions.csv"))

cat("Done. Figures and predictions saved to:", OUT_DIR, "/\n")
