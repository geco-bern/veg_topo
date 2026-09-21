# ============================================================================
# 1. LIBRARIES AND CONFIGURATION
# ============================================================================

library(dplyr)
library(tidyr)
library(caret)
library(recipes)
library(ggplot2)
library(arrow)
library(ranger)
library(patchwork)

source(here::here("R/config.R"))

# Set parallel cores
if (hostname == "dash") n_cores = 32 else n_cores = 80
message("→ using ", n_cores, " workers\n")

# Check if output directory exists, create if not
if (!dir.exists(rf_models_dir)) dir.create(rf_models_dir, recursive = TRUE)

# ============================================================================
# 2. LOAD DATA
# ============================================================================

# Load data
dataset <- open_dataset(rf_sample_data_tiles_dir, format = "parquet")
combined_data <- dataset |> collect() |> tidyr::drop_na() |> filter(lat <= 60)
message("→ Loaded ", nrow(combined_data), " samples\n")
head(combined_data)
# ============================================================================
# 3. FIXED PARAMETERS FOR ALL MODELS (CRITICAL FOR VARIANCE PARTITIONING)
# ============================================================================

FIXED_PARAMS <- list(
  mtry = 2,               # Fixed across all models (balanced randomness)
  min.node.size = 10,     # Fixed
  num.trees = 500,        # Fixed
  sample.fraction = 0.8,  # Fixed
  replace = FALSE,
  importance = "permutation"
)

message("=== FIXED PARAMETERS (identical for all models) ===")
message(sprintf("  mtry = %d", FIXED_PARAMS$mtry))
message(sprintf("  min.node.size = %d", FIXED_PARAMS$min.node.size))
message(sprintf("  num.trees = %d", FIXED_PARAMS$num.trees))
message(sprintf("  sample.fraction = %.1f\n", FIXED_PARAMS$sample.fraction))

# ============================================================================
# 4. FUNCTION TO TRAIN MODEL WITH FIXED PARAMETERS
# ============================================================================

train_rf_model <- function(formula, data, model_name, n_cores) {

  message("\n", paste(rep("=", 50), collapse = ""))
  message("Training: ", model_name)
  message(paste(rep("=", 50), collapse = ""))

  # Print sample info
  message(sprintf("  Training samples: %d", nrow(data)))
  message(sprintf("  Features: %d", length(all.vars(formula)) - 1))  # 减去因变量
  message(sprintf("  Model formula: %s", deparse(formula)))

  # Create recipe
  pp <- recipe(formula, data = data) |>
    step_center(all_numeric(), -all_outcomes()) |>
    step_scale(all_numeric(), -all_outcomes())

  # Training control
  train_control <- trainControl(
    method = "cv",
    number = 5,
    savePredictions = "final",
    verboseIter = FALSE,
    allowParallel = FALSE
  )

  # Train model with FIXED parameters
  set.seed(1982)
  start_time <- Sys.time()

  model <- train(
    pp,
    data = data,
    method = "ranger",
    trControl = train_control,
    tuneGrid = data.frame(
      mtry = FIXED_PARAMS$mtry,
      splitrule = "variance",
      min.node.size = FIXED_PARAMS$min.node.size
    ),
    metric = "RMSE",
    replace = FIXED_PARAMS$replace,
    sample.fraction = FIXED_PARAMS$sample.fraction,
    num.trees = FIXED_PARAMS$num.trees,
    importance = FIXED_PARAMS$importance,
    num.threads = n_cores,
    seed = 1982
  )

  end_time <- Sys.time()
  training_time <- difftime(end_time, start_time, units = "mins")

  # Extract results
  cv_rmse <- model$results$RMSE
  cv_rsq <- model$results$Rsquared
  oob_rmse <- sqrt(model$finalModel$prediction.error)

  # Variable importance
  importance <- data.frame(
    Variable = names(model$finalModel$variable.importance),
    Importance = model$finalModel$variable.importance
  ) |> arrange(desc(Importance))

  message(sprintf("  ✓ Completed in %.2f minutes", training_time))
  message(sprintf("  ✓ CV R² = %.4f", cv_rsq))
  message(sprintf("  ✓ CV RMSE = %.4f", cv_rmse))
  message(sprintf("  ✓ OOB RMSE = %.4f", oob_rmse))

  return(list(
    model = model,
    name = model_name,
    cv_rmse = cv_rmse,
    cv_rsq = cv_rsq,
    oob_rmse = oob_rmse,
    importance = importance,
    predictions = model$pred,
    training_time = training_time
  ))
}

# ============================================================================
# 5. TRAIN THREE MODELS WITH IDENTICAL PARAMETERS
# ============================================================================

cat("\n", paste(rep("=", 60), collapse = ""))
cat("\nTRAINING THREE MODELS FOR VARIANCE PARTITIONING\n")
cat(paste(rep("=", 60), collapse = ""), "\n")


# Model A: Full model (Climate + Topography)
gc()
full <- train_rf_model(
  vegh ~ mat + map + srad + twi + rin + elv,
  combined_data,
  "Full Model (Climate + Topo)",
  n_cores
)
# Save individual models
saveRDS(full$model, file.path(rf_models_dir, paste0("rf_full_model_cv.rds")))


# Model B: Climate only
gc()
climate <- train_rf_model(
  vegh ~ mat + map + srad,
  combined_data,
  "Climate Only",
  n_cores
)
saveRDS(climate$model, file.path(rf_models_dir, paste0("rf_climate_model_cv.rds")))

# Model C: Topography only
gc()
topo <- train_rf_model(
  vegh ~ twi + rin + elv,
  combined_data,
  "Topography Only",
  n_cores
)
saveRDS(topo$model, file.path(rf_models_dir, paste0("rf_topo_model_cv.rds")))

# ============================================================================
# 6. CROSS-VALIDATION RESULTS FOR EACH MODEL
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""))
cat("\nCROSS-VALIDATION RESULTS\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Full model results
cat("\n", paste(rep("-", 40), collapse = ""))
cat("\nFull Model (Climate + Topography)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for(i in 1:5) {
  fold_preds <- full$predictions[full$predictions$Resample == paste0("Fold", i), ]
  rmse <- sqrt(mean((fold_preds$pred - fold_preds$obs)^2))
  rsq <- cor(fold_preds$pred, fold_preds$obs)^2
  cat(sprintf("Fold %d: RMSE = %.4f, R² = %.4f\n", i, rmse, rsq))
}
cat(sprintf("Average: RMSE = %.4f, R² = %.4f\n", full$cv_rmse, full$cv_rsq))

# Climate only results
cat("\n", paste(rep("-", 40), collapse = ""))
cat("\nClimate Only Model\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for(i in 1:5) {
  fold_preds <- climate$predictions[climate$predictions$Resample == paste0("Fold", i), ]
  rmse <- sqrt(mean((fold_preds$pred - fold_preds$obs)^2))
  rsq <- cor(fold_preds$pred, fold_preds$obs)^2
  cat(sprintf("Fold %d: RMSE = %.4f, R² = %.4f\n", i, rmse, rsq))
}
cat(sprintf("Average: RMSE = %.4f, R² = %.4f\n", climate$cv_rmse, climate$cv_rsq))

# Topography only results
cat("\n", paste(rep("-", 40), collapse = ""))
cat("\nTopography Only Model\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
for(i in 1:5) {
  fold_preds <- topo$predictions[topo$predictions$Resample == paste0("Fold", i), ]
  rmse <- sqrt(mean((fold_preds$pred - fold_preds$obs)^2))
  rsq <- cor(fold_preds$pred, fold_preds$obs)^2
  cat(sprintf("Fold %d: RMSE = %.4f, R² = %.4f\n", i, rmse, rsq))
}
cat(sprintf("Average: RMSE = %.4f, R² = %.4f\n", topo$cv_rmse, topo$cv_rsq))

# ============================================================================
# 7. VARIABLE IMPORTANCE FOR EACH MODEL
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""))
cat("\nVARIABLE IMPORTANCE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("\n--- Full Model ---\n")
print(full$importance)

cat("\n--- Climate Only Model ---\n")
print(climate$importance)

cat("\n--- Topography Only Model ---\n")
print(topo$importance)

# ============================================================================
# 8. VARIANCE PARTITIONING
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""))
cat("\nVARIANCE PARTITIONING RESULTS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

R2_full <- full$cv_rsq
R2_climate <- climate$cv_rsq
R2_topo <- topo$cv_rsq

# Calculate variance components
unique_climate <- R2_full - R2_topo
unique_topo <- R2_full - R2_climate
joint <- R2_climate + R2_topo - R2_full
unexplained <- 1 - R2_full

cat("Model R² values (all with identical parameters):\n")
cat(sprintf("  Full Model (Climate + Topo):  %.4f\n", R2_full))
cat(sprintf("  Climate Only:                  %.4f\n", R2_climate))
cat(sprintf("  Topography Only:               %.4f\n\n", R2_topo))

cat("Variance Components:\n")
cat(sprintf("  Unique Climate contribution:   %.4f (%.1f%% of explained)\n",
            unique_climate, unique_climate / R2_full * 100))
cat(sprintf("  Unique Topography contribution: %.4f (%.1f%% of explained)\n",
            unique_topo, unique_topo / R2_full * 100))
cat(sprintf("  Joint (Climate + Topo):        %.4f (%.1f%% of explained)\n",
            joint, joint / R2_full * 100))
cat(sprintf("  Total explained:               %.4f (%.1f%% of total)\n",
            R2_full, R2_full * 100))
cat(sprintf("  Unexplained:                   %.4f (%.1f%% of total)\n\n",
            unexplained, unexplained * 100))

# ============================================================================
# 9. VISUALIZATION: PREDICTED VS OBSERVED (THREE MODELS)
# ============================================================================

# Prepare data for plotting
full_pred <- full$predictions[, c("pred", "obs", "Resample")]
full_pred$Model <- "Full Model"

climate_pred <- climate$predictions[, c("pred", "obs", "Resample")]
climate_pred$Model <- "Climate Only"

topo_pred <- topo$predictions[, c("pred", "obs", "Resample")]
topo_pred$Model <- "Topography Only"

all_pred <- rbind(full_pred, climate_pred, topo_pred)

# Calculate R² for each model
model_rsq <- data.frame(
  Model = c("Full Model", "Climate Only", "Topography Only"),
  R2 = c(R2_full, R2_climate, R2_topo),
  RMSE = c(full$cv_rmse, climate$cv_rmse, topo$cv_rmse)
)

# Plot 1: Predicted vs observed by model
p1 <- ggplot(all_pred, aes(x = obs, y = pred, color = Resample)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1) +
  facet_wrap(~Model, nrow = 1) +
  labs(
    title = "Predicted vs Observed (5-fold CV)",
    x = "Observed", y = "Predicted"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# print(p1)

# Plot 2: R² comparison
p2 <- ggplot(model_rsq, aes(x = Model, y = R2, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.4f", R2)), vjust = -0.5, size = 4) +
  labs(
    title = "Model Performance Comparison",
    subtitle = "All models use identical parameters (mtry=2, num.trees=500)",
    y = expression(R^2)
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(0, 1)

# print(p2)

# Plot 3: Variance partitioning pie/bar chart
variance_data <- data.frame(
  Component = c("Unique Climate", "Unique Topography", "Joint", "Unexplained"),
  Value = c(unique_climate, unique_topo, joint, unexplained)
)

p3 <- ggplot(variance_data, aes(x = "", y = Value, fill = Component)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", Value * 100)),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  labs(
    title = "Variance Partitioning",
    subtitle = paste0("Total R² = ", round(R2_full, 4))
  ) +
  theme_void() +
  theme(legend.position = "bottom")

# print(p3)

# Plot 4: Variable importance comparison
imp_combined <- rbind(
  cbind(full$importance, Model = "Full"),
  cbind(climate$importance, Model = "Climate Only"),
  cbind(topo$importance, Model = "Topography Only")
)

p4 <- ggplot(imp_combined, aes(x = reorder(Variable, Importance),
                               y = Importance, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_wrap(~Model, scales = "free_y", nrow = 1) +
  labs(
    title = "Variable Importance by Model",
    x = "Variables", y = "Importance (Permutation)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

# print(p4)

# Save plots
ggsave(here::here("data/figures/7_02_pred_vs_obs_three_models.png"),
       p1, width = 15, height = 5, dpi = 300)
ggsave(here::here("data/figures/7_02_r2_comparison.png"),
       p2, width = 6, height = 6, dpi = 300)
ggsave(here::here("data/figures/7_02_variance_partitioning.png"),
       p3, width = 8, height = 6, dpi = 300)
ggsave(here::here("data/figures/7_02_variable_importance_comparison.png"),
       p4, width = 15, height = 5, dpi = 300)

# ============================================================================
# 10. SAVE RESULTS
# ============================================================================

# Save variance partitioning results
variance_summary <- data.frame(
  Component = c("Unique Climate", "Unique Topography", "Joint", "Total Explained", "Unexplained"),
  R2 = c(unique_climate, unique_topo, joint, R2_full, unexplained),
  Percentage_of_Total = c(
    unique_climate, unique_topo, joint, R2_full, unexplained
  ) * 100,
  Percentage_of_Explained = c(
    unique_climate / R2_full * 100,
    unique_topo / R2_full * 100,
    joint / R2_full * 100,
    100,
    NA
  )
)

model_comparison <- data.frame(
  Model = c("Full", "Climate Only", "Topography Only"),
  CV_RMSE = c(full$cv_rmse, climate$cv_rmse, topo$cv_rmse),
  CV_R2 = c(full$cv_rsq, climate$cv_rsq, topo$cv_rsq),
  OOB_RMSE = c(full$oob_rmse, climate$oob_rmse, topo$oob_rmse),
  Training_Time_Min = c(full$training_time, climate$training_time, topo$training_time)
)


# Save all results
saveRDS(model_comparison, file.path(rf_models_dir, "model_comparison_cv.rds"))
saveRDS(variance_summary, file.path(rf_models_dir, "variance_summary_cv.rds"))
saveRDS(full$importance, file.path(rf_models_dir, "importance_full_cv.rds"))
saveRDS(climate$importance, file.path(rf_models_dir, "importance_climate_cv.rds"))
saveRDS(topo$importance, file.path(rf_models_dir, "importance_topo_cv.rds"))

# Also save as CSV for easy viewing
write.csv(model_comparison, file.path(rf_models_dir, "model_comparison_cv.csv"), row.names = FALSE)
write.csv(variance_summary, file.path(rf_models_dir, "variance_summary_cv.csv"), row.names = FALSE)

# Save combined importance for visualization
importance_combined <- bind_rows(
  full$importance %>% mutate(Model = "Full"),
  climate$importance %>% mutate(Model = "Climate Only"),
  topo$importance %>% mutate(Model = "Topography Only")
)
saveRDS(importance_combined, file.path(rf_models_dir, "importance_combined_cv.rds"))
# ============================================================================
# 11. TRAINING SUMMARY
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""))
cat("\nTRAINING SUMMARY\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Total training time:\n")
cat(sprintf("  Full Model:        %.2f minutes\n", full$training_time))
cat(sprintf("  Climate Only:      %.2f minutes\n", climate$training_time))
cat(sprintf("  Topography Only:   %.2f minutes\n", topo$training_time))
cat(sprintf("  Total:             %.2f minutes\n",
            full$training_time + climate$training_time + topo$training_time))

cat("\nParameters used (identical for all models):\n")
cat(sprintf("  mtry = %d\n", FIXED_PARAMS$mtry))
cat(sprintf("  num.trees = %d\n", FIXED_PARAMS$num.trees))
cat(sprintf("  min.node.size = %d\n", FIXED_PARAMS$min.node.size))
cat(sprintf("  sample.fraction = %.1f\n", FIXED_PARAMS$sample.fraction))


cat("\n✅ COMPLETED!\n")
