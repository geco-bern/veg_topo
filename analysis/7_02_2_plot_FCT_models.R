# ============================================================================
# VISUALIZATION: THREE MODELS WITH R² AND VARIABLE IMPORTANCE
# ============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(here)
library(khroma)

source(here::here("R/config.R"))

# ============================================================================
# 1. LOAD SAVED MODEL RESULTS
# ============================================================================

# Load model comparison and importance data
model_comparison <- readRDS(file.path(rf_models_dir, "model_comparison.rds"))
importance_full <- readRDS(file.path(rf_models_dir, "importance_full.rds"))
importance_climate <- readRDS(file.path(rf_models_dir, "importance_climate.rds"))
importance_topo <- readRDS(file.path(rf_models_dir, "importance_topo.rds"))

# ============================================================================
# 2. PREPARE DATA - COMBINE R² WITH VARIABLE IMPORTANCE
# ============================================================================

# Calculate percentages for each model's importance and sort by importance
importance_full_pct <- importance_full %>%
  mutate(
    Percentage = Importance / sum(Importance) * 100,
    Model = "Full Model"
  ) %>%
  arrange(desc(Percentage))

importance_climate_pct <- importance_climate %>%
  mutate(
    Percentage = Importance / sum(Importance) * 100,
    Model = "Climate Only"
  ) %>%
  arrange(desc(Percentage))

importance_topo_pct <- importance_topo %>%
  mutate(
    Percentage = Importance / sum(Importance) * 100,
    Model = "Topography Only"
  ) %>%
  arrange(desc(Percentage))

# IMPORTANT: Create a factor with order from largest to smallest
# This controls the stacking order in ggplot
full_order <- importance_full_pct$Variable  # Order from largest to smallest

# For each model, ensure Variable is a factor with the correct order
importance_full_pct$Variable <- factor(importance_full_pct$Variable, levels = full_order)
importance_climate_pct$Variable <- factor(importance_climate_pct$Variable, levels = full_order)
importance_topo_pct$Variable <- factor(importance_topo_pct$Variable, levels = full_order)

# Combine all data
all_importance <- bind_rows(
  importance_full_pct,
  importance_climate_pct,
  importance_topo_pct
)

# Add R² and RMSE values to the data
r2_values <- model_comparison %>%
  mutate(
    Model = case_when(
      Model == "Full" ~ "Full Model",
      Model == "Climate Only" ~ "Climate Only",
      Model == "Topography Only" ~ "Topography Only"
    )
  ) %>%
  select(Model, R2 = Training_R2, RMSE = Training_RMSE)

# Merge importance with R²
all_importance <- all_importance %>%
  left_join(r2_values, by = "Model") %>%
  mutate(
    R2_contribution = (Percentage / 100) * R2
  )

# ============================================================================
# 3. SET MODEL ORDER (Top to Bottom: Full, Climate Only, Topography Only)
# ============================================================================

# Define model order - this controls the order on Y-axis (top to bottom)
model_order <- c("Topography Only", "Climate Only", "Full Model")

# Convert Model to factor with specified order
all_importance$Model <- factor(all_importance$Model, levels = model_order)
r2_values$Model <- factor(r2_values$Model, levels = model_order)

# ============================================================================
# 4. CREATE HORIZONTAL STACKED BAR CHART
# ============================================================================

# Define color palette for variables (ordered by importance in Full Model)

# 创建调色板并提取 6 个颜色
batlow_pal <- color("batlow")
batlow_colors <- batlow_pal(6)

# 应用到变量
variable_colors <- c(
  "map" = batlow_colors[1],
  "mat" = batlow_colors[2],
  "srad" = batlow_colors[3],
  "elv" = batlow_colors[4],
  "rin" = batlow_colors[5],
  "twi" = batlow_colors[6]
)
# Make sure colors match the factor levels
variable_colors <- variable_colors[names(variable_colors) %in% levels(all_importance$Variable)]

# Create horizontal stacked bar plot
# For horizontal bars with coord_flip():
# - position_stack() with reverse = FALSE: bottom-to-top in legend = left-to-right in bars
# - We want largest importance on the left, so we reverse the factor levels for stacking
text_size = 7
p_horizontal <- ggplot(all_importance, aes(x = Model, y = R2_contribution, fill = Variable)) +

  # Stacked bars - reverse = TRUE makes the first factor level (largest importance) go to bottom/left
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE),
           width = 0.6, color = "white", linewidth = 0.1) +
  # # Add percentage labels inside each segment
  # geom_text(
  #   aes(
  #     label = ifelse(Percentage > 3, sprintf("%.1f%%", Percentage), ""),
  #     group = Variable
  #   ),
  #   position = position_stack(vjust = 0.5, reverse = TRUE),
  #   size = 1.5,
  #   color = "black",
  #   hjust = 0.5,
  #   vjust = 0.5,        # 垂直对齐，配合旋转
  #   angle = 90,         # 旋转90度
  #   lineheight = 0.9
  # ) +

  # Add total R² and RMSE at the end of each bar (to the right)
  geom_text(
    data = r2_values,
    aes(x = Model, y = R2 + 0.02,
        label = sprintf("R² = %.4f\nRMSE = %.2f m", R2, RMSE)),
    inherit.aes = FALSE,
    size = 2,
    color = "black",
    hjust = 0,
    lineheight = 0.9
  ) +

  # Use manual colors
  scale_fill_manual(
    values = variable_colors,
    name = "",
    labels = c(
      "map" = "MAP",
      "mat" = "MAT",
      "srad" = "SRAD",
      "elv" = "Elv",
      "rin" = "Rin",
      "twi" = "TWI"
    ),
    # This ensures the legend shows the same order
    breaks = full_order
  ) +

  # Flip coordinates to make horizontal bars
  coord_flip() +

  scale_x_discrete(
    labels = c(
      "Full Model" = "F",
      "Climate Only" = "C",
      "Topography Only" = "T"
    )
  ) +

  # Labels and theme
  labs(
    # title = "Model Performance with Variable Importance Contributions",
    # subtitle = expression("Bar length" ~ R^2 ~ "; Colors from left to right show decreasing importance"),
    x = NULL,
    y = expression(R^2)
  ) +
  theme_bw(base_size = text_size) +
  theme(
    plot.title = element_blank(),

    # 坐标轴标签
    axis.text.y = element_text(size = text_size, face = "bold"),
    axis.text.x = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),

    # 图例设置 - 大幅减小
    legend.position = "right",
    legend.direction = "vertical",  # 垂直排列
    legend.box = "vertical",        # 垂直盒子
    legend.text = element_text(size = text_size),  # 从 text_size 减小到 0.8倍
    legend.title = element_text(size = text_size),
    legend.key.size = unit(0.1, "in"),  # 减小图例色块大小
    legend.key.width = unit(0.1, "in"),  # 减小图例色块宽度
    legend.key.height = unit(0.22, "in"),  # 减小图例色块高度
    legend.spacing = unit(0.1, "in"),  # 减小图例间距
    legend.margin = margin(0, 0, 0, 0, "pt"),  # 减小图例边距

    # 面板网格
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线


    # 整体边距
    plot.margin = margin(1, 1, 1, 1, "mm")  # 减小图形边距
  ) +

  scale_y_continuous(
    limits = c(0, max(r2_values$R2) * 1.30),  # 从 1.35 减小到 1.30
    expand = expansion(mult = c(0, 0.02))  # 减小扩展
  )

# ============================================================================
# 5. SAVE THE PLOT
# ============================================================================

ggsave(
  here::here("data/figures/7_02_model_r2_horizontal_with_importance.png"),
  p_horizontal,
  width = 3.5,  # Slightly wider to accommodate RMSE text
  height = 1.8,
  dpi =300
)

# ============================================================================
# 6. SUMMARY OUTPUT
# ============================================================================

cat("\n", paste(rep("=", 60), collapse = ""))
cat("\nVARIABLE ORDER (from largest to smallest importance)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("Full Model order (left to right in bars):\n")
print(full_order)

cat("\nModel order (top to bottom):\n")
print(model_order)

cat("\nIn horizontal bars: leftmost (closest to Y-axis) = largest importance\n")
cat("Models arranged from top to bottom: Full → Climate Only → Topography Only\n\n")

cat("Model Performance Metrics:\n")
print(r2_values)

cat("\n✅ Done! Saved to: data/figures/7_02_model_r2_horizontal_with_importance.png\n")
