# =============================================================================
# Setup - Load required packages
# =============================================================================
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tidyterra)
library(scico)
library(rnaturalearth)
library(sf)
library(khroma)

source(here::here("R/config.R"))
source(here::here("R/plot_cor_pval.R"))

# =============================================================================
# 1. Data preparation
# =============================================================================

# 1.1 Load coastline data
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# 1.2 Read Rin-vegetation height correlation raster
r_H_Rin <- rast(r_H_R_5km_path)

# 1.3 Convert to data frame and clean
r_H_Rin_clean <- as.data.frame(r_H_Rin, xy = FALSE) |>
  drop_na()
colnames(r_H_Rin_clean)[1] <- "correlation"

# 1.4 Create grouping variable (positive / negative correlation)
r_H_Rin_clean <- r_H_Rin_clean |>
  mutate(group = case_when(
    correlation > 0 ~ "Positive",
    correlation < 0 ~ "Negative"
  ))

# 1.5 Calculate percentages of positive and negative correlations
percent_positive <- sum(r_H_Rin_clean$correlation > 0) / nrow(r_H_Rin_clean) * 100
percent_negative <- sum(r_H_Rin_clean$correlation < 0) / nrow(r_H_Rin_clean) * 100

# 1.6 Pre-calculate maximum bin count for annotation positioning
bin_count <- max(table(cut(r_H_Rin_clean$correlation, breaks = seq(-1, 1, by = 0.1))))

# 1.7 Extract colors from vik palette (red for positive, green for negative)
palette_colors <- scico(256, palette = "vik")
positive_color <- palette_colors[50]  # 红色
negative_color <- palette_colors[200]   # 蓝色

# =============================================================================
# 2. Create correlation histogram (p1)
# =============================================================================
p1 <- ggplot(r_H_Rin_clean, aes(x = correlation, fill = group)) +
  # Histogram: binwidth 0.05, black border, semi-transparent
  geom_histogram(binwidth = 0.05, color = "black", linewidth = 0.2, alpha = 0.5, boundary = 0) +
  # Manual fill colors
  scale_fill_manual(
    values = c("Positive" = positive_color, "Negative" = negative_color),
    name = "Correlation Type"
  ) +
  # Dashed vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  # x-axis limits
  xlim(-1, 1) +
  # Axis labels
  labs(
    title = NULL,
    x = "r(H,Rᵢₙ)",
    y = expression("Count (×10"^4*")")
  ) +
  # Divide y-axis labels by 10000
  scale_y_continuous(
    labels = function(x) round(x / 10000, 1)
  ) +
  # Classic theme with custom styling
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
    # All font sizes set to 7
    axis.title.x = element_text(hjust = 0.5, margin = margin(t = 1), size = 7),
    axis.title.y = element_text(hjust = 0.5, angle = 90, margin = margin(r = 1), size = 7),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_line(linewidth = 0.3),
    axis.ticks.y = element_line(linewidth = 0.3),
    axis.line = element_line(color = "black", linewidth = 0.2),
    panel.background = element_blank()
  ) +
  # Percentage annotations
  annotate("text", x = 0.5, y = 400000,
           label = sprintf("%.1f%%", percent_positive),
           color = palette_colors[50], size = 2.5,
           fontface = "bold",
           vjust = -0.5) +
  annotate("text", x = -0.5, y = 400000,
           label = sprintf("%.1f%%", percent_negative),
           color = palette_colors[200], size = 2.5,
           fontface = "bold",
           vjust = -0.5)

# =============================================================================
# 3. Create global correlation map (p_cor)
# =============================================================================

# 3.1 Aggregate raster for coarser resolution
r_H_Rin <- aggregate(r_H_Rin, fact = c(2, 2))

# 3.2 Build map
p_cor <- ggplot2::ggplot() +
  # Spatial raster layer
  tidyterra::geom_spatraster(data = r_H_Rin, maxcell = Inf) +
  # Color scale using bam palette
  scico::scale_fill_scico(
    palette = "vik",
    direction = -1,
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5),
    midpoint = 0,
    na.value = NA
  ) +
  # Labels
  ggplot2::labs(
    title = NULL,
    fill = expression(r(H, R["in"])),
  ) +
  # x-axis (longitude)
  ggplot2::scale_x_continuous(
    breaks = seq(from = -180, to = 180, by = 30),
    limits = c(-180, 180),
    expand = expansion(mult = 0.0001)
  ) +
  # y-axis (latitude)
  ggplot2::scale_y_continuous(
    breaks = seq(from = -60, to = 90, by = 30),
    limits = c(-60, 90),
    expand = expansion(mult = 0.0001)
  ) +
  # Theme
  ggplot2::theme_bw(base_size = 7) +
  # Colorbar guide
  guides(fill = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = grid::unit(0.1, "in"),
    barheight = grid::unit(2.6, "in")
  )) +
  # Coastline overlay
  geom_sf(data = coast,
          colour = 'black',
          linewidth = 0.1) +
  # Theme customization
  ggplot2::theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -6),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.5),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.spacing = unit(0, "pt"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )


# =============================================================================
# 4. Add tags and combine plots
# =============================================================================

# 4.1 Add tag "b)" to histogram
p1 <- p1 +
  labs(tag = "b)") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.1, 1.1),
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_blank()
  )

# 4.2 Add tag "a)" to map
p_cor <- p_cor +
  labs(tag = "a)") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_blank()
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_blank()
  )

# 4.3 Combine: map + inset histogram
p_cor_comb <- p_cor +
  inset_element(
    p1 +
      theme(
        plot.background = element_rect(fill = NA, color = NA),
        plot.margin = margin(0, 0, 0, 0)
      ),
    left = 0.022,
    bottom = 0.05,
    right = 0.27,
    top = 0.45,
    align_to = "full"
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_blank()
  )


ggsave(
  filename = file.path(project_root, "data/figures/3_01_r_H_Rin_5km_map_0.1d.png"),
  plot = p_cor_comb,
  width = 7,
  height = 3,
  dpi = 600,
  units = "in"
)


# =============================================================================
# 5. Plot p-value map (p_pval)
# =============================================================================

# 5.1 Aggregate p-value raster
p_H_Rin <- rast(pval_r_H_R_5km_path) |>
  aggregate(fact = c(2, 2))

# 5.2 Create p-value map using custom function
p_pval <- plot_cor_pval(
  input = p_H_Rin,
  extent = ext_global,
  title_text = NULL,
  fill_text = "r(H,Rᵢₙ)",
  text_size = 7,
  x_step = 30,
  y_step = 30
) +
  # Colorbar guide
  guides(fill = guide_colorbar(
    title.position = "top",
    barwidth = grid::unit(0.1, "in"),
    barheight = grid::unit(2.6, "in")
  )) +
  # Coastline overlay
  geom_sf(data = coast,
          colour = 'black',
          linewidth = 0.1) +
  # Theme customization
  ggplot2::theme(
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -6),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5, vjust = 0.5),
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, vjust = 0.5),
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    panel.spacing = unit(0, "pt"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank()
  )

ggsave(
  filename = file.path(project_root, "data/figures/3_01_r_H_Rin_5km_pval_0.1d.png"),
  plot = p_pval,
  width = 7,
  height = 3,
  dpi = 600,
  units = "in"
)
