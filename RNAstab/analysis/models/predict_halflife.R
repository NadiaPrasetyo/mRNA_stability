# predict_halflife.R
#
# Apply a trained half-life model to a feature dataframe.
# The model .rds must have been produced by train_halflife_model.R.
#
# Usage in a script:
#   source("predict_halflife.R")
#   bovine_df <- readRDS("bovine_features.rds")
#   preds     <- predict_halflife("models/halflife_lgbm_v1.rds", bovine_df)
#   write_csv(preds, "bovine_halflife_predictions.csv")
#
# Required packages:
#   install.packages(c("tidyverse", "tidymodels", "bonsai", "lightgbm"))

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(bonsai)
})

predict_halflife <- function(model_path,
                             feature_df,
                             return_id_cols = TRUE) {

  if (!file.exists(model_path)) {
    stop("Model file not found: ", model_path)
  }

  m <- readRDS(model_path)

  # ----- Validate that all required feature columns are present -----
  required <- m$feature_cols
  missing  <- setdiff(required, names(feature_df))

  if (length(missing) > 0) {
    stop(
      "Input dataframe is missing ", length(missing),
      " required feature column(s):\n  ",
      paste(head(missing, 30), collapse = "\n  "),
      if (length(missing) > 30) "\n  (and more...)" else "",
      "\n\nThe model was trained with these features and cannot",
      " predict without them. Compute them with the upstream feature pipeline,",
      " or fill with NA if the column genuinely doesn't apply."
    )
  }

  extra <- setdiff(
    names(feature_df),
    c(required, m$meta_cols, m$target_col)
  )
  if (length(extra) > 0) {
    message("Note: ", length(extra),
            " input columns not used by the model (ignored).")
  }

  # ----- Predict -----
  preds <- predict(m$workflow, feature_df) %>%
    rename(predicted_halflife = .pred)

  # ----- Attach IDs back if available -----
  id_present <- intersect(m$meta_cols, names(feature_df))
  if (return_id_cols && length(id_present) > 0) {
    preds <- bind_cols(feature_df[id_present], preds)
  }

  # ----- Report ------
  rsq <- m$test_metrics %>%
    filter(.metric == "rsq") %>%
    pull(.estimate)

  message(sprintf(
    "Predicted %d transcripts. Model trained on %d %s rows; held-out R2 = %.3f, Spearman = %.3f.",
    nrow(preds), m$n_train, m$trained_on, rsq, m$test_spearman
  ))

  preds
}


# ---------- Example invocation (uncomment and edit) ----------
# bovine_df <- readRDS("bovine_features.rds")
# preds     <- predict_halflife("models/halflife_lgbm_v1.rds", bovine_df)
# write_csv(preds, "bovine_halflife_predictions.csv")
