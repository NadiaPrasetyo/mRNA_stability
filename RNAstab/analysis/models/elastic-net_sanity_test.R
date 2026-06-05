# baseline_glmnet.R
#
# Linear (elastic net) baseline for the half-life model. Same data, same split,
# same target — only the model differs. Compare test R² against the LightGBM
# model's test R² to decide whether the tree complexity is justified.

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(finetune)
  library(glmnet)        # install.packages("glmnet") if needed
})

set.seed(42)             # MUST match training script for split reproducibility

INPUT_PATH  <- "data/cache/human_dataset_v4.rds"
OUTPUT_PATH <- "models/halflife_glmnet_baseline_v4.rds"
N_CORES     <- max(1, parallel::detectCores() - 1)

META_COLS   <- c("species", "transcript_id", "gene_id", "gene_name")
TARGET_COL  <- "halflife"

# Same drop logic as the LightGBM training script — keep these in sync, or
# better, factor them into a shared config later.
DROP_PATTERNS <- c("^rnafold_", "^rnalfold_", "^rnaup_",
                   "^mfe_(expected|delta)_", "^(shape|keth)_", "^length_")
EXTRA_DROP <- c("orfexondensity", "expression", "translation_efficiency",
                "junctions_count_mrna", "cds_length_codons_cds", "aa_total_cds",
                "n_stops_cds", "n_codons_scored_cds", "uorf_count_mrna",
                "utr5_length", "junctions_count_cds", "junctions_count_3utr",
                "junctions_count_5utr", "n_exons", "exon_count_internal_mrna",
                "cai")

df <- readRDS(INPUT_PATH)
drop_cols <- df %>%
  select(matches(paste(DROP_PATTERNS, collapse = "|"))) %>% names()

df_model <- df %>%
  select(-all_of(c(drop_cols, EXTRA_DROP))) %>%
  filter(!is.na(.data[[TARGET_COL]]))

feature_cols <- setdiff(names(df_model), c(META_COLS, TARGET_COL))

split    <- initial_split(df_model, prop = 0.8)
train_df <- training(split)
test_df  <- testing(split)

# Linear-specific recipe: NA imputation + standardisation. Keep step_corr too
# so the comparison is on equal footing with your latest LightGBM recipe.
glmnet_recipe <- recipe(train_df) %>%
  update_role(all_of(feature_cols), new_role = "predictor") %>%
  update_role(all_of(TARGET_COL),   new_role = "outcome") %>%
  update_role(any_of(META_COLS),    new_role = "id") %>%
  step_zv(all_predictors()) %>%
  step_impute_median(all_predictors()) %>%
  step_corr(all_predictors(), threshold = 0.85) %>%
  step_normalize(all_predictors())

glmnet_spec <- linear_reg(
  penalty = tune(),    # lambda (overall regularisation strength)
  mixture = tune()     # alpha: 0 = ridge, 1 = lasso, in between = elastic net
) %>%
  set_engine("glmnet") %>%
  set_mode("regression")

glmnet_wf <- workflow() %>%
  add_recipe(glmnet_recipe) %>%
  add_model(glmnet_spec)

glmnet_grid <- grid_space_filling(
  penalty(range = c(-4, 0)),    # log10 scale: 1e-4 to 1
  mixture(range = c(0, 1)),
  size = 25
)

folds <- vfold_cv(train_df, v = 5)

cat("--- Tuning elastic net baseline ---\n")
tune_res <- tune_race_anova(
  glmnet_wf,
  resamples = folds,
  grid      = glmnet_grid,
  metrics   = metric_set(rmse, rsq),
  control   = control_race(verbose_elim = TRUE)
)

best_params <- select_best(tune_res, metric = "rmse")
final_fit   <- finalize_workflow(glmnet_wf, best_params) %>% fit(train_df)

# After final_fit, before saveRDS:
coefs <- final_fit %>%
  extract_fit_engine() %>%
  coef(s = best_params$penalty) %>%
  as.matrix() %>%
  as_tibble(rownames = "feature") %>%
  rename(coefficient = `s1`) %>%
  filter(coefficient != 0, feature != "(Intercept)") %>%
  arrange(desc(abs(coefficient)))

cat(sprintf("\nLasso retained %d of %d features.\n",
            nrow(coefs), length(feature_cols)))
cat("Top 15 by |coefficient|:\n")
print(head(coefs, 15))

test_preds   <- augment(final_fit, test_df)
test_metrics <- test_preds %>% metrics(truth = !!sym(TARGET_COL), estimate = .pred)
spearman     <- cor(test_preds[[TARGET_COL]], test_preds$.pred,
                    method = "spearman", use = "complete.obs")

cat("\n--- Elastic net test performance ---\n")
print(test_metrics)
cat(sprintf("Spearman: %.3f\n", spearman))
cat("\nBest hyperparameters:\n"); print(best_params)

saveRDS(
  list(workflow = final_fit, best_params = best_params,
       test_metrics = test_metrics, test_spearman = spearman,
       feature_cols = feature_cols, meta_cols = META_COLS,
       target_col = TARGET_COL),
  OUTPUT_PATH
)
