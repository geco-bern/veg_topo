library(terra)
library(tidyverse)
library(patchwork)
source(here::here("R/config.R"))
source(here::here("R/mosaic_tiles.R"))

# 定义处理函数（不返回数字，只保存图形）
process_and_plot <- function(input_dir, pattern_lidar, pattern_lang, output_name) {

  message("⭐️⭐️⭐️ Processing: ", output_name, " ⭐️⭐️⭐️")

  # 读取栅格
  r_lidar <- mosaic_tiles(input_dir = input_dir, output_file = NULL, pattern = pattern_lidar)
  r_lang  <- mosaic_tiles(input_dir = input_dir, output_file = NULL, pattern = pattern_lang)

  # 堆叠并转换
  stacked <- c(r_lidar, r_lang)
  names(stacked) <- c("lidar", "lang")

  df <- as.data.frame(stacked, xy = TRUE, na.rm = TRUE)
  df$cor_diff <- df$lang - df$lidar

  # 添加统计信息到图形标题（可选）
  cor_coef <- cor(df$lang, df$lidar, use = "complete.obs")
  mean_diff <- mean(df$cor_diff, na.rm = TRUE)

  # 定义统一的主题模板（增大字体）
  base_theme <- theme_minimal() +
    theme(
      # 全局字体大小
      text = element_text(size = 14),
      # 坐标轴标题字体
      axis.title = element_text(size = 16, face = "bold"),
      # 坐标轴刻度标签字体
      axis.text = element_text(size = 12),
      # 图例文字字体
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      # 标题字体
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      # 分面标题（如果有）
      strip.text = element_text(size = 14, face = "bold")
    )

  # 创建三个图形
  p1 <- ggplot(df, aes(x = lang, y = lidar)) +
    geom_hex(bins = 30) +
    scale_fill_viridis_c(name = "Count", trans = "log") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", linewidth = 1) +
    labs(title = paste(output_name, ": Correlation Comparison"),
         subtitle = paste("R =", round(cor_coef, 3)),
         x = "Lang Correlation", y = "Lidar Correlation") +
    base_theme +
    theme(legend.position = "none")

  p2 <- ggplot(df, aes(x = cor_diff)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = mean_diff, color = "darkblue", linetype = "dotted", linewidth = 1) +
    labs(title = paste(output_name, ": Difference Distribution"),
         subtitle = paste("Mean diff =", round(mean_diff, 4)),
         x = "Difference (Lang - Lidar)") +
    base_theme

  p3 <- df %>%
    pivot_longer(cols = c(lang, lidar), names_to = "variable", values_to = "correlation") %>%
    ggplot(aes(x = variable, y = correlation, fill = variable)) +
    geom_boxplot() +
    scale_fill_manual(values = c("lang" = "lightgreen", "lidar" = "lightcoral"),
                      labels = c("Lang", "Lidar")) +
    labs(title = paste(output_name, ": Boxplot Comparison"),
         y = "Correlation", x = "") +
    base_theme +
    theme(legend.position = "none")

  # 4. 密度图（两个图层叠加）
  p4 <- df %>%
    pivot_longer(cols = c(lang, lidar), names_to = "variable", values_to = "correlation") %>%
    ggplot(aes(x = correlation, fill = variable, color = variable)) +
    geom_density(alpha = 0.5, linewidth = 0.8) +
    scale_fill_manual(values = c("lang" = "lightgreen", "lidar" = "lightcoral"),
                      labels = c("Lang", "Lidar")) +
    scale_color_manual(values = c("lang" = "darkgreen", "lidar" = "darkred"),
                       labels = c("Lang", "Lidar")) +
    labs(title = "Density Distribution",
         x = "Correlation", y = "Density") +
    base_theme +
    theme(legend.position = "right",
          legend.title = element_blank())

  # 组合图形 2x2 布局
  combined_plot <- (p1 + p2) / (p3 + p4) +
    plot_annotation(title = paste(output_name, "Analysis"),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))

  # 确保输出目录存在
  output_dir <- here::here("data/figures")
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("📁 Created directory: ", output_dir)
  }

  # 保存组合图
  output_file <- file.path(output_dir, paste0("6_01_validate_H_", output_name, "_all.png"))
  ggsave(output_file, combined_plot, width = 12, height = 8, dpi = 300)
  message("✅ Figure saved: ", output_file)
  message("   📈 Correlation (R): ", round(cor_coef, 3))
  message("   📊 Mean difference: ", round(mean_diff, 5))

  # 不返回任何值
  invisible(NULL)
}

# 定义参数列表
params <- list(
  TWI_500m = list(
    pattern_lidar = "*_r_hlidar_twi.tif",
    pattern_lang  = "*_r_hlang_twi.tif",
    name = "TWI_500m"
  ),
  TWI_5km = list(
    pattern_lidar = "*_r_hlidar_twi_5km.tif",
    pattern_lang  = "*_r_hlang_twi_5km.tif",
    name = "TWI_5km"
  ),
  Rin_500m = list(
    pattern_lidar = "*_r_hlidar_rin.tif",
    pattern_lang  = "*_r_hlang_rin.tif",
    name = "Rin_500m"
  ),
  Rin_5km = list(
    pattern_lidar = "*_r_hlidar_rin_5km.tif",
    pattern_lang  = "*_r_hlang_rin_5km.tif",
    name = "Rin_5km"
  )
)

# 循环处理并保存图形
for (param in params) {
  process_and_plot(
    input_dir = h_validation_dir,
    pattern_lidar = param$pattern_lidar,
    pattern_lang = param$pattern_lang,
    output_name = param$name
  )
}

message("\n🎉 All plots have been saved to data/figures/ directory!")
