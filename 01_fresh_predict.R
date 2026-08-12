library(dplyr)
library(caret)

# 0. Data

data_spp <- read.csv("all_trait_data.csv")
df_fresh_dry <- read.csv("data_fresh_dry.csv")

all_vars <- names(df_fresh_dry)

base_traits <- unique(gsub("_fresh|_dry", "", all_vars))

# Keep only complete pairs

base_traits <- base_traits[sapply(base_traits, function(v) {
  all(c(paste0(v, "_fresh"), paste0(v, "_dry")) %in% all_vars)
})]

# 1. TRAINING (LM)

metrics_list <- list()
saved_models <- list()

set.seed(123)

for (trait in base_traits) {
  
  fresh_col <- paste0(trait, "_fresh")
  dry_col   <- paste0(trait, "_dry")
  
  model_data <- df_fresh_dry %>%
    dplyr::select(all_of(fresh_col), all_of(dry_col)) %>%
    na.omit()
  
  if (nrow(model_data) >= 20) {
    
    message("Processing trait: ", trait)
    
    # Linear model
    formula <- as.formula(paste(fresh_col, "~", dry_col))
    lm_model <- lm(formula, data = model_data)
    
    # Holdout
    train_id <- createDataPartition(
      model_data[[fresh_col]],
      p = 0.8,
      list = FALSE
    )
    
    pred <- predict(lm_model, newdata = model_data[-train_id, ])
    obs  <- model_data[[fresh_col]][-train_id]
    
    rmse <- sqrt(mean((obs - pred)^2))
    r2   <- cor(obs, pred)^2
    mean_val <- mean(obs)
    range_obs <- max(obs) - min(obs)
    
    saved_models[[trait]] <- list(
      model = lm_model,
      predictor = dry_col
    )
    
    metrics_list[[trait]] <- data.frame(
      trait = trait,
      sample_size = nrow(model_data),
      rmse = rmse,
      r2 = r2,
      mean_fresh_val = mean_val,
      range = range_obs,
      nrmse = rmse / range_obs
    )
  }
}

summary_performance <- bind_rows(metrics_list)

print(summary_performance)

write.csv(
  summary_performance,
  file = "fresh_prediction_lm_metrics.csv",
  row.names = FALSE
)

# 2. Convert mm to cm

cols_to_convert <- c(
  "petal_size_mm", "flower_diam_mm", "flower_length_mm",
  "hipant_diam_mm", "hipant_length_mm", "fruit_length_mean_mm",
  "fruit_diameter_mean_mm", "seed_length_mean_mm", "seed_diameter_mean_mm"
)

data_spp <- data_spp %>%
  mutate(across(
    all_of(cols_to_convert),
    ~ . / 10,
    .names = "{.col}_cm"
  ))

# 3. Mapping

column_map <- list(
  "flower_petal_size"      = "petal_size_mm_cm",
  "flower_diameter_petal"  = "flower_diam_mm_cm",
  "flower_length"          = "flower_length_mm_cm",
  "flower_hipant_diameter" = "hipant_diam_mm_cm",
  "flower_hipant_length"   = "hipant_length_mm_cm",
  "fruit_length"            = "fruit_length_mean_mm_cm",
  "fruit_diameter"          = "fruit_diameter_mean_mm_cm",
  "seed_length"             = "seed_length_mean_mm_cm",
  "seed_diameter"           = "seed_diameter_mean_mm_cm",
  
  # Already in cm
  
  "leaf_area"              = "leaf_area_cm2",
  "leaf_perimeter"         = "leaf_perim_cm",
  "leaf_length"            = "leaf_length_cm",
  "leaf_width"             = "leaf_width_cm",
  "petiole_length"         = "petiole_length_cm",
  "petiole_width_base"     = "petiole_width_at_base_cm",
  "petiole_width_middle"   = "petiole_width_at_50pct_cm",
  "petiole_width_top"      = "petiole_width_at_top_cm",
  "leaf_length_per_width"   = "leaf_length_per_width",
  "leaf_circulariry_index" = "circularity_index",
  "leaf_lma_area"           = "lma_pred"
)

# 4. Prediction (LM)

data_spp_final <- data_spp

for (trait in names(saved_models)) {
  
  model_obj <- saved_models[[trait]]
  target_col <- column_map[[trait]]
  pred_name  <- paste0(trait, "_fresh_pred")
  
  if (!is.null(target_col) && target_col %in% names(data_spp_final)) {
    
    message("Predicting: ", trait)
    
    data_spp_final[[pred_name]] <- NA
    
    valid_idx <- which(!is.na(data_spp_final[[target_col]]))
    
    if (length(valid_idx) > 0) {
      
      newdata <- data_spp_final[
        valid_idx,
        target_col,
        drop = FALSE
      ]
      
      # Rename to the model predictor name
      colnames(newdata) <- model_obj$predictor
      
      data_spp_final[[pred_name]][valid_idx] <-
        predict(model_obj$model, newdata = newdata)
    }
    
  } else {
    warning("Skipping trait: ", trait)
  }
}

# 5. Replace values

data_spp_fresh <- data_spp

for (trait in names(saved_models)) {
  
  data_col <- column_map[[trait]]
  pred_col <- paste0(trait, "_fresh_pred")
  
  if (!(pred_col %in% names(data_spp_final))) next
  
  fresh_values <- data_spp_final[[pred_col]]
  
  # Convert back to mm if necessary
  
  if (grepl("*mm*", data_col)) {
    fresh_values <- fresh_values * 10
  }
  
  original_name <- gsub("_cm$", "", data_col)
  
  if (!(original_name %in% names(data_spp_fresh))) {
    original_name <- data_col
  }
  
  new_name <- paste0(original_name, "_fresh")
  
  data_spp_fresh[[new_name]] <- fresh_values
  
  message("Created: ", new_name)
}

# 6. Cleaning

data_spp_fresh <- data_spp_fresh %>%
  select(-ends_with("_cm"))

# 7. Output

write.csv(
  data_spp_fresh,
  file = "all_trait_data_fresh_predicted_LM.csv",
  row.names = FALSE
)

message(
  "Done! Final dimensions: ",
  nrow(data_spp_fresh),
  " x ",
  ncol(data_spp_fresh)
)