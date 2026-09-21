# ============================================================================
# 1. LIBRARIES AND CONFIGURATION
# ============================================================================

library(dplyr)
library(tidyr)
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
combined_data <- dataset |>
  collect() |>
  tidyr::drop_na() |>
  filter(lat <= 60)
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
# 4. FUNCTION TO TRAIN MODEL WITH FIXED PARAMETERS (USING RANGER)
# ============================================================================

train_rf_model_ranger <- function(formula, data, model_name, n_cores) {

  message("\n", paste(rep("=", 50), collapse = ""))
  message("Training: ", model_name)
  message(paste(rep("=", 50), collapse = ""))

  # Print sample info
  message(sprintf("  Training samples: %d", nrow(data)))
  message(sprintf("  Features: %d", length(all.vars(formula)) - 1))
  message(sprintf("  Model formula: %s", deparse(formula)))

  set.seed(1982)
  start_time <- Sys.time()

  # Core training - one line with ranger
  model <- ranger(
    formula,
    data = data,
    num.trees = FIXED_PARAMS$num.trees,
    mtry = FIXED_PARAMS$mtry,
    min.node.size = FIXED_PARAMS$min.node.size,
    importance = FIXED_PARAMS$importance,
    replace = FIXED_PARAMS$replace,
    sample.fraction = FIXED_PARAMS$sample.fraction,
    num.threads = n_cores
  )

  end_time <- Sys.time()
  training_time <- difftime(end_time, start_time, units = "mins")

  # Extract outcome name from formula
  outcome_name <- as.character(formula[[2]])
  observed <- data[[outcome_name]]

  # Calculate training RMSE
  train_rmse <- sqrt(mean((model$predictions - observed)^2))

  message(sprintf("  ✓ Completed in %.2f minutes", training_time))
  message(sprintf("  ✓ Training R² = %.4f", model$r.squared))
  message(sprintf("  ✓ Training RMSE = %.4f", train_rmse))
  message(sprintf("  ✓ OOB RMSE = %.4f", sqrt(model$prediction.error)))

  # Return results
  return(list(
    model = model,
    name = model_name,
    oob_rmse = sqrt(model$prediction.error),
    train_rsq = model$r.squared,
    train_rmse = train_rmse,
    predictions = model$predictions,
    importance = model$variable.importance,
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
full <- train_rf_model_ranger(
  vegh ~ mat + map + srad + twi + rin + elv,
  combined_data,
  "Full Model (Climate + Topo)",
  n_cores
)
saveRDS(full$model, file.path(rf_models_dir, "rf_full_model.rds"))

# Model B: Climate only
gc()
climate <- train_rf_model_ranger(
  vegh ~ mat + map + srad,
  combined_data,
  "Climate Only",
  n_cores
)
saveRDS(climate$model, file.path(rf_models_dir, "rf_climate_model.rds"))

# Model C: Topography only
gc()
topo <- train_rf_model_ranger(
  vegh ~ twi + rin + elv,
  combined_data,
  "Topography Only",
  n_cores
)
saveRDS(topo$model, file.path(rf_models_dir, "rf_topo_model.rds"))

# ============================================================================
# 6. MODEL PERFORMANCE (USING OOB ERROR)
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""))
cat("\nMODEL PERFORMANCE (OOB ERROR)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Full model results
cat("\n", paste(rep("-", 40), collapse = ""))
cat("\nFull Model (Climate + Topography)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat(sprintf("OOB RMSE = %.4f\n", full$oob_rmse))
cat(sprintf("Training R² = %.4f\n", full$train_rsq))
cat(sprintf("Training RMSE = %.4f\n", full$train_rmse))

# Climate only results
cat("\n", paste(rep("-", 40), collapse = ""))
cat("\nClimate Only Model\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat(sprintf("OOB RMSE = %.4f\n", climate$oob_rmse))
cat(sprintf("Training R² = %.4f\n", climate$train_rsq))
cat(sprintf("Training RMSE = %.4f\n", climate$train_rmse))

# Topography only results
cat("\n", paste(rep("-", 40), collapse = ""))
cat("\nTopography Only Model\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat(sprintf("OOB RMSE = %.4f\n", topo$oob_rmse))
cat(sprintf("Training R² = %.4f\n", topo$train_rsq))
cat(sprintf("Training RMSE = %.4f\n", topo$train_rmse))

# ============================================================================
# 7. VARIABLE IMPORTANCE FOR EACH MODEL
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""))
cat("\nVARIABLE IMPORTANCE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("\n--- Full Model ---\n")
importance_full <- data.frame(
  Variable = names(full$importance),
  Importance = as.numeric(full$importance)
) |> arrange(desc(Importance))
print(importance_full)

cat("\n--- Climate Only Model ---\n")
importance_climate <- data.frame(
  Variable = names(climate$importance),
  Importance = as.numeric(climate$importance)
) |> arrange(desc(Importance))
print(importance_climate)

cat("\n--- Topography Only Model ---\n")
importance_topo <- data.frame(
  Variable = names(topo$importance),
  Importance = as.numeric(topo$importance)
) |> arrange(desc(Importance))
print(importance_topo)

# ============================================================================
# 8. VARIANCE PARTITIONING (USING TRAINING R²)
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""))
cat("\nVARIANCE PARTITIONING RESULTS\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

R2_full <- full$train_rsq
R2_climate <- climate$train_rsq
R2_topo <- topo$train_rsq

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
full_pred <- data.frame(
  pred = full$predictions,
  obs = combined_data$vegh,
  Model = "Full Model"
)

climate_pred <- data.frame(
  pred = climate$predictions,
  obs = combined_data$vegh,
  Model = "Climate Only"
)

topo_pred <- data.frame(
  pred = topo$predictions,
  obs = combined_data$vegh,
  Model = "Topography Only"
)

all_pred <- rbind(full_pred, climate_pred, topo_pred)

# Calculate R² and RMSE for each model for plotting
model_rsq <- data.frame(
  Model = c("Full Model", "Climate Only", "Topography Only"),
  R2 = c(full$train_rsq, climate$train_rsq, topo$train_rsq),
  RMSE = c(full$train_rmse, climate$train_rmse, topo$train_rmse)
)

# Plot 1: Predicted vs observed by model
p1 <- ggplot(all_pred, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.4, size = 1, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1, color = "red") +
  facet_wrap(~Model, nrow = 1) +
  labs(
    title = "Predicted vs Observed (Training Set)",
    x = "Observed", y = "Predicted"
  ) +
  theme_bw()

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

# Plot 4: Variable importance comparison
# Combine importance data from all models
full_imp_df <- data.frame(
  Variable = names(full$importance),
  Importance = as.numeric(full$importance),
  Model = "Full"
)

climate_imp_df <- data.frame(
  Variable = names(climate$importance),
  Importance = as.numeric(climate$importance),
  Model = "Climate Only"
)

topo_imp_df <- data.frame(
  Variable = names(topo$importance),
  Importance = as.numeric(topo$importance),
  Model = "Topography Only"
)

imp_combined <- rbind(full_imp_df, climate_imp_df, topo_imp_df)

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
  Training_RMSE = c(full$train_rmse, climate$train_rmse, topo$train_rmse),
  Training_R2 = c(full$train_rsq, climate$train_rsq, topo$train_rsq),
  OOB_RMSE = c(full$oob_rmse, climate$oob_rmse, topo$oob_rmse),
  Training_Time_Min = c(full$training_time, climate$training_time, topo$training_time)
)

# Save results
saveRDS(variance_summary, file.path(rf_models_dir, "variance_partitioning.rds"))
saveRDS(model_comparison, file.path(rf_models_dir, "model_comparison.rds"))
saveRDS(importance_full, file.path(rf_models_dir, "importance_full.rds"))
saveRDS(importance_climate, file.path(rf_models_dir, "importance_climate.rds"))
saveRDS(importance_topo, file.path(rf_models_dir, "importance_topo.rds"))

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
cat(sprintf("  replace = %s\n", FIXED_PARAMS$replace))
cat(sprintf("  importance = %s\n", FIXED_PARAMS$importance))

cat("\n✅ COMPLETED!\n")
