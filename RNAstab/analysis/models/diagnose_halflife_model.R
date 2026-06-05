# diagnose_halflife_model.R
#
# Diagnostic plots for the half-life model produced by train_halflife_model.R.
#
# Reproduces the 80/20 split (set.seed(42)) used in training, refits on the
# training portion using the saved best hyperparameters, predicts on the
# held-out test set, and produces:
#
#   Dashboard (2x3, saved as diagnostic_dashboard.{png,pdf}):
#     1. Predicted vs observed (hexbin)
#     2. Residuals vs predicted
#     3. Residual distribution
#     4. Calibration by predicted decile
#     5. Observed half-life within each predicted decile (violin + box)
#     6. Top-K recovery curve (filtering performance)
#
#   Standalone:
#     7. Top 25 features by gain
#
# Also saves the test-set predictions to figures/test_predictions.csv.
#
# Required packages:
#   install.packages(c("tidyverse", "tidymodels", "bonsai", "lightgbm",
#                      "patchwork", "hexbin"))
#
# Usage:
#   Edit the paths below, then: Rscript diagnose_halflife_model.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(bonsai)
  library(patchwork)
})

set.seed(42)  # must match training script

# ----------------------------- Configuration --------------------------------

INPUT_PATH  <- "data/cache/human_dataset_v5.rds"            # <-- training dataframe
MODEL_PATH  <- "models/halflife_lgbm_v5.rds"        # <-- saved model
FIGURES_DIR <- "figures"
N_CORES     <- max(1, parallel::detectCores() - 1)

dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)


# ----------------------------- 1. Load and reproduce the split --------------

m  <- readRDS(MODEL_PATH)
df <- readRDS(INPUT_PATH)

target_col   <- m$target_col
feature_cols <- m$feature_cols

# Use the saved feature set rather than re-deriving drops here. Filtering down
# also keeps initial_split() reproducible across schema changes upstream.
df_model <- df %>%
  select(any_of(m$meta_cols), all_of(feature_cols), all_of(target_col)) %>%
  filter(!is.na(.data[[target_col]]))

split    <- initial_split(df_model, prop = 0.8)
train_df <- training(split)
test_df  <- testing(split)


# ----------------------------- 2. Refit on train, predict on test -----------

cat("Refitting on training set with saved hyperparameters...\n")

hl_recipe <- recipe(train_df) %>%
  update_role(all_of(feature_cols), new_role = "predictor") %>%
  update_role(all_of(target_col),   new_role = "outcome") %>%
  update_role(any_of(m$meta_cols),  new_role = "id") %>%
  step_zv(all_predictors())

bp <- m$best_params

lgbm_spec <- boost_tree(
  trees          = 1000,
  learn_rate     = bp$learn_rate,
  tree_depth     = bp$tree_depth,
  min_n          = bp$min_n,
  loss_reduction = bp$loss_reduction,
  sample_size    = bp$sample_size,
  mtry           = bp$mtry
) %>%
  set_engine("lightgbm", num_threads = N_CORES, counts = FALSE) %>%
  set_mode("regression")

train_fit <- workflow() %>%
  add_recipe(hl_recipe) %>%
  add_model(lgbm_spec) %>%
  fit(train_df)

test_pred <- augment(train_fit, test_df) %>%
  filter(!is.na(.pred), !is.na(.data[[target_col]])) %>%
  mutate(residual = .pred - .data[[target_col]])

train_pred <- augment(train_fit, train_df) %>%
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

s_test  <- stats(test_pred,  target_col)
s_train <- stats(train_pred, target_col)

cat(sprintf("\nTrain: n=%d, R^2=%.3f, Spearman=%.3f, RMSE=%.3f\n",
            s_train$n, s_train$rsq, s_train$rho, s_train$rmse))
cat(sprintf("Test:  n=%d, R^2=%.3f, Spearman=%.3f, RMSE=%.3f\n\n",
            s_test$n,  s_test$rsq,  s_test$rho,  s_test$rmse))

saved_rmse <- m$test_metrics %>%
  filter(.metric == "rmse") %>% pull(.estimate)

if (abs(s_test$rmse - saved_rmse) > 0.01) {
  warning(sprintf(
    "Refit RMSE (%.3f) differs from saved test RMSE (%.3f) — pipeline drift?",
    s_test$rmse, saved_rmse
  ))
}

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
# Switch geom_hex() to geom_point(alpha = 0.05) if you'd rather see raw points.

p_truth <- ggplot(test_pred, aes(x = .data[[target_col]], y = .pred)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  geom_abline(slope = 1, intercept = 0, colour = "white",
              linetype = "dashed", linewidth = 0.6) +
  geom_smooth(method = "lm", colour = "firebrick", se = FALSE,
              linewidth = 0.6, formula = y ~ x) +
  labs(title    = "Predicted vs observed (held-out)",
       subtitle = stat_label(s_test),
       x        = "Observed (PC1 half-life)",
       y        = "Predicted") +
  theme_diag()


# ----------------------------- 6. Plot 2: Residuals vs predicted ------------

p_resid <- ggplot(test_pred, aes(x = .pred, y = residual)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "count") +
  geom_hline(yintercept = 0, colour = "white",
             linetype = "dashed", linewidth = 0.6) +
  geom_smooth(method = "loess", colour = "firebrick", se = FALSE,
              linewidth = 0.6, formula = y ~ x) +
  labs(title    = "Residuals vs predicted",
       subtitle = "Loess curve should hover around zero",
       x        = "Predicted",
       y        = "Residual (predicted \u2212 observed)") +
  theme_diag()


# ----------------------------- 7. Plot 3: Residual distribution -------------

resid_combined <- bind_rows(
  train_pred %>% mutate(split = "Train"),
  test_pred  %>% mutate(split = "Test")
) %>%
  mutate(split = factor(split, levels = c("Train", "Test")))

p_rdist <- ggplot(resid_combined, aes(x = residual, fill = split, colour = split)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed",
             linewidth = 0.4) +
  scale_fill_manual(values   = c(Train = "grey55", Test = "steelblue"), name = NULL) +
  scale_colour_manual(values = c(Train = "grey55", Test = "steelblue"), name = NULL) +
  labs(
    title    = "Residual distribution: train vs test",
    subtitle = sprintf(
      "Train: mean = %.3f, sd = %.3f  |  Test: mean = %.3f, sd = %.3f",
      mean(train_pred$residual), sd(train_pred$residual),
      mean(test_pred$residual),  sd(test_pred$residual)
    ),
    x = "Residual (predicted \u2212 observed)",
    y = "Density"
  ) +
  theme_diag() +
  theme(legend.position = "top")

# ----------------------------- 8. Plot 4: Calibration -----------------------

cal_data <- test_pred %>%
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


# ----------------------------- 9. Plot 5: Observed within predicted decile --
# Direct view of bin separability — overlap = poor ranking quality

p_perfbin <- test_pred %>%
  mutate(bin = ntile(.pred, 10)) %>%
  ggplot(aes(x = factor(bin), y = .data[[target_col]])) +
  geom_violin(fill = "steelblue", alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  labs(title    = "Observed half-life within predicted decile",
       subtitle = "Less overlap between adjacent bins = better ranking",
       x        = "Predicted decile (1 = lowest, 10 = highest)",
       y        = "Observed (PC1 half-life)") +
  theme_diag()


# ----------------------------- 10. Plot 6: Top-K recovery -------------------
# If you take the top K% by predicted half-life, what % of the true top K%
# do you capture?  Direct measure of filtering quality.

ks <- seq(0.01, 1, by = 0.01)

recovery <- map_dfr(ks, function(k) {
  cutoff_pred  <- quantile(test_pred$.pred,            1 - k)
  cutoff_truth <- quantile(test_pred[[target_col]],    1 - k)
  predicted_top <- test_pred$.pred         >= cutoff_pred
  true_top      <- test_pred[[target_col]] >= cutoff_truth
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


# Trim to 5th-95th percentile of observed — precision/recall get noisy at the
# tails where the positive class is tiny.
cutoffs <- seq(
  quantile(test_pred[[target_col]], 0.05, na.rm = TRUE),
  quantile(test_pred[[target_col]], 0.95, na.rm = TRUE),
  length.out = 60
)

threshold_data <- map_dfr(cutoffs, function(c) {
  predicted_pos <- test_pred$.pred         >= c
  true_pos      <- test_pred[[target_col]] >= c
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
    title    = "Filtering quality at absolute cutoffs",
    subtitle = "Same cutoff applied to predicted and observed half-life score",
    x        = "Half-life score cutoff",
    y        = "%"
  ) +
  theme_diag() +
  theme(legend.position = "bottom")

# ----------------------------- 11. Combine into dashboard -------------------

dashboard <- (p_truth | p_resid | p_rdist) /
             (p_cal   | p_perfbin | p_topk) +
  plot_annotation(
    title    = "Half-life model diagnostic dashboard",
    subtitle = sprintf(
      "Held-out (n = %d) vs training (n = %d) | R\u00B2 test = %.3f, train = %.3f",
      s_test$n, s_train$n, s_test$rsq, s_train$rsq
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(colour = "grey40")
    )
  )

ggsave(file.path(FIGURES_DIR, "diagnostic_dashboard.png"),
       dashboard, width = 16, height = 10, dpi = 200)
ggsave(file.path(FIGURES_DIR, "diagnostic_dashboard.pdf"),
       dashboard, width = 16, height = 10)

filter_combined <- (p_topk | p_thresh) +
  plot_annotation(
    title    = "Filtering quality — relative vs absolute cutoffs",
    subtitle = "Same model, two framings: 'top K%' (rank-based) vs 'score \u2265 c' (threshold-based)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey40")
    )
  )

ggsave(file.path(FIGURES_DIR, "filtering_quality.png"),
       filter_combined, width = 13, height = 5.5, dpi = 200)
ggsave(file.path(FIGURES_DIR, "filtering_quality.pdf"),
       filter_combined, width = 13, height = 5.5)

# ----------------------------- 12. Plot 7: Feature importance ---------------

if (!is.null(m$importance)) {
  imp_df <- m$importance %>%
    as_tibble() %>%
    slice_head(n = 25) %>%
    mutate(Feature = fct_reorder(Feature, Gain))

  p_imp <- ggplot(imp_df, aes(x = Gain, y = Feature)) +
    geom_col(fill = "steelblue") +
    labs(title    = "Top 25 features by gain",
         subtitle = "From the LightGBM final model fit on full data",
         x        = "Gain (relative)",
         y        = NULL) +
    theme_diag() +
    theme(panel.grid.major.y = element_blank())

  ggsave(file.path(FIGURES_DIR, "feature_importance.png"),
         p_imp, width = 8, height = 8, dpi = 200)
  ggsave(file.path(FIGURES_DIR, "feature_importance.pdf"),
         p_imp, width = 8, height = 8)
} else {
  warning("No importance table found in model .rds — skipping plot 7.")
}


# ----------------------------- 13. Save predictions for reuse ---------------

test_pred %>%
  select(any_of(m$meta_cols), all_of(target_col), .pred, residual) %>%
  write_csv(file.path(FIGURES_DIR, "test_predictions.csv"))

cat("Done. Figures and predictions saved to:", FIGURES_DIR, "/\n")
