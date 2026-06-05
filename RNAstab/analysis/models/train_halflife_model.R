# train_halflife_model.R
#
# Trains a LightGBM model to predict mRNA half-life from sequence-derived
# features. Designed for cross-species deployment (currently: bovine).
#
# Target: PC1 of consensus half-life from Agarwal & Kelley 2022.
# Training data: human only (one row per gene, MANE-select transcript).
#
# Required packages:
#   install.packages(c("tidyverse", "tidymodels", "bonsai", "finetune", "lightgbm"))
#
# Usage:
#   1. Edit INPUT_PATH and OUTPUT_PATH below
#   2. Rscript train_halflife_model.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(bonsai)        # LightGBM engine for parsnip
  library(finetune)      # tune_race_anova for fast tuning
})

set.seed(42)

# ----------------------------- Configuration --------------------------------

INPUT_PATH  <- "data/cache/human_dataset_v5.rds"
OUTPUT_PATH <- "models/halflife_lgbm_v5.rds"
N_CORES     <- max(1, parallel::detectCores() - 1)

# Columns that are IDs / metadata / target — never used as predictors.
META_COLS  <- c("species", "transcript_id", "gene_id", "gene_name")
TARGET_COL <- "halflife"

# Patterns to drop:
#   - Vienna features (rnafold, rnalfold, rnaup, mfe_*) are slow on bovine
#   - Probing features (shape, keth) are human-only
DROP_PATTERNS <- c(
  "^rnafold_",
  "^rnalfold_",
  "^rnaup_",
  "^mfe_(expected|delta)_",
  "^(shape|keth)_",
  "^length_"
)

# Anything else to manually exclude (column names, not patterns)
EXTRA_DROP <- c("orfexondensity", "expression", "translation_efficiency", "junctions_count_mrna", 
                "cds_length_codons_cds", "aa_total_cds", "n_stops_cds", "n_codons_scored_cds",
                "uorf_count_mrna", "utr5_length", "junctions_count_cds","junctions_count_3utr",
                "junctions_count_5utr", "n_exons", "exon_count_internal_mrna", "cai", "stop_dist_closest_downstream",
                "stop_dist_last_downstream", "max_classical_uorf_codons", "dist_last_uorf_stop_to_main_atg",
                "start_dist_closest_upstream", "dist_cap_to_first_uatg_mrna", "internal_exon_sd", "internal_exon_mean",
                "internal_exon_median", "intron_sd")


# ----------------------------- 1. Load and shape ----------------------------

df <- readRDS(INPUT_PATH)

drop_cols <- df %>%
  select(matches(paste(DROP_PATTERNS, collapse = "|"))) %>%
  names()

df_model <- df %>%
  select(-all_of(c(drop_cols, EXTRA_DROP))) %>%
  filter(!is.na(.data[[TARGET_COL]]))

feature_cols <- setdiff(names(df_model), c(META_COLS, TARGET_COL))

cat("\n--- Setup ---\n")
cat("Rows after dropping missing target:", nrow(df_model), "\n")
cat("Predictor features:                ", length(feature_cols), "\n")
cat("Dropped (Vienna / probing):        ", length(drop_cols), "\n\n")


# ----------------------------- 2. Train / test split ------------------------
# Random 80/20. If you have chromosome info, switch to grouping by chromosome
# via group_initial_split() to be defensible against paralog leakage.

split    <- initial_split(df_model, prop = 0.8)
train_df <- training(split)
test_df  <- testing(split)


# ----------------------------- 3. Recipe ------------------------------------
# LightGBM tolerates NAs natively (good — bovine inputs will have them too),
# and is scale-invariant, so the recipe is minimal: just role assignment plus
# a zero-variance filter.

hl_recipe <- recipe(train_df) %>%
  update_role(all_of(feature_cols), new_role = "predictor") %>%
  update_role(all_of(TARGET_COL),   new_role = "outcome") %>%
  update_role(any_of(META_COLS),    new_role = "id") %>%
  step_nzv(all_predictors()) %>%
  step_corr(all_predictors(), threshold = 0.75)


# ----------------------------- 4. Model spec --------------------------------
lgbm_spec <- boost_tree(
  trees          = 500,                    # was 1000 — pair with the tightened ranges below
  learn_rate     = tune(),
  tree_depth     = tune(),
  min_n          = tune(),
  loss_reduction = tune(),
  sample_size    = tune(),
  mtry           = tune()
) %>%
  set_engine("lightgbm",
             num_threads = N_CORES,
             counts      = FALSE,
             lambda_l2   = tune("lambda_l2")) %>%   # L2 weight penalty on leaf values
  set_mode("regression")


# ----------------------------- 5. Tuning ------------------------------------

hl_wf <- workflow() %>%
  add_recipe(hl_recipe) %>%
  add_model(lgbm_spec)



hl_params <- hl_wf %>%
  extract_parameter_set_dials() %>%
  update(
    mtry           = mtry_prop(c(0.3, 1.0)),
    learn_rate     = learn_rate(c(-3, -1.5)),       # was c(-3, -1) — cap a touch lower
    tree_depth     = tree_depth(c(3, 8)),           # was c(3, 12) — cap depth
    min_n          = min_n(c(20, 100)),             # was c(5, 40) — raise floor and ceiling
    loss_reduction = loss_reduction(c(-2, 1)),
    sample_size    = sample_prop(c(0.5, 1.0)),
    lambda_l2      = dials::new_quant_param(
      type      = "double",
      range     = c(-4, 1),        # log10 scale: 1e-4 to 10
      inclusive = c(TRUE, TRUE),
      trans     = scales::log10_trans(),
      label     = c(lambda_l2 = "L2 regularisation")
    )
  )

folds <- vfold_cv(train_df, v = 5)
lgbm_grid <- grid_space_filling(hl_params, size = 30)  # 30 instead of 25 to absorb the extra dim


cat("--- Tuning ---\n")
tune_res <- tune_race_anova(
  hl_wf,
  resamples = folds,
  grid      = lgbm_grid,
  metrics   = metric_set(rmse, rsq),
  control   = control_race(verbose_elim = TRUE, save_pred = FALSE)
)


# ----------------------------- 6. Finalise ----------------------------------

best_params <- select_best(tune_res, metric = "rmse")
final_wf    <- finalize_workflow(hl_wf, best_params)

cat("\n--- Best hyperparameters ---\n")
print(best_params)

cat("\n--- Fitting on training set for evaluation ---\n")
final_fit <- fit(final_wf, train_df)


# ----------------------------- 7. Held-out evaluation -----------------------

test_preds <- augment(final_fit, test_df)

test_metrics <- test_preds %>%
  metrics(truth = !!sym(TARGET_COL), estimate = .pred)

spearman <- cor(test_preds[[TARGET_COL]], test_preds$.pred,
                method = "spearman", use = "complete.obs")

cat("\n--- Held-out test performance ---\n")
print(test_metrics)
cat(sprintf("Spearman rho: %.3f\n", spearman))


# ----------------------------- 8. Refit on all data and save ----------------

cat("\n--- Refitting on full dataset for deployment ---\n")
final_full_fit <- fit(final_wf, df_model)


# ----------------------------- 9. Feature importance (sanity check) ---------

imp <- final_full_fit %>%
  extract_fit_engine() %>%
  lightgbm::lgb.importance() %>%
  as_tibble()

cat("\n--- Top 20 features by gain ---\n")
print(imp %>% slice_head(n = 20))


# ----------------------------- 10. Save -------------------------------------

dir.create(dirname(OUTPUT_PATH), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  list(
    workflow      = final_full_fit,
    feature_cols  = feature_cols,
    meta_cols     = META_COLS,
    target_col    = TARGET_COL,
    best_params   = best_params,
    test_metrics  = test_metrics,
    test_spearman = spearman,
    importance    = imp,
    trained_on    = "human (Agarwal & Kelley PC1)",
    n_train       = nrow(df_model),
    fitted_at     = Sys.time(),
    notes         = paste(
      "LightGBM regressor for mRNA half-life.",
      "Designed for cross-species deployment.",
      "Inference dataframe must contain the same feature column names."
    )
  ),
  OUTPUT_PATH
)

cat("\nSaved model to:", OUTPUT_PATH, "\n")
