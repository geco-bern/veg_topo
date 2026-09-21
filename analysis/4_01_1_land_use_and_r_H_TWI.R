# ==============================================================================
# Script: ~/veg_topo/analysis/4_01_1_land_use_and_r_H_TWI.R
#
# Purpose:
#   Combined analysis of TWI-VEGH correlation patterns:
#     1. Mask correlation data into natural and human-used lands
#     2. Create scatter plots across land-use intensity bins
#     3. Compare correlation distributions between land types (natural and human-used lands)
#     4. Combine all visualizations in a comprehensive layout
#
# Dependencies:
#   - terra, dplyr, ggplot2, viridis, patchwork, tidyterra, scico, rnaturalearth, sf
#   - config.R or config_ubelix.R (for file paths)
#   - R/plot_hex_scatter.R, R/plot_cor_twi_vegh.R (custom plotting functions)
# ==============================================================================

# ------------------------- 1. Setup Environment --------------------------------
library(terra)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(tidyterra)
library(scico)
library(rnaturalearth)
library(sf)
library(grid) # add figure tag

source(here::here("R/config.R"))
source(here::here("R/plot_hex_scatter.R"))
source(here::here("R/plot_cor_twi_vegh.R"))

# ------------------------- 2. Load and Prepare Data ----------------------------
# Coastline (for map overlay)
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# Load rasters
fused_5km_r <- rast(fused_5km_file)
mi_5km_r <- rast(mi_5km_file) * 0.0001
cor_twi_vegh_r <- rast(cor_twi_vegh_mosaic_file)

# Apply masks to correlation data for distribution comparison
cor_used <- mask(cor_twi_vegh_r, fused_5km_r >= 0.05, maskvalues = FALSE) # correlation on human-used land
cor_natural <- mask(cor_twi_vegh_r, fused_5km_r < 0.05, maskvalues = FALSE) # correlation on natural land

# ------------------------- 3. Plot Global Correlation Maps (单独保存) ----------
## 3.1. Used land (fused ≥ 0.05)
p_used <- plot_cor_twi_vegh(
  input = cor_used,
  extent = ext_global,
  title = "Pearson’s r (H~TWI) on used land (fused ≥ 0.05)",
  text_size = 12,
  x_step = 30, y_step = 30
) +
  guides(fill = guide_colorbar(
    title.position = "left",
    barwidth = grid::unit(0.1, "in"),
    barheight = grid::unit(5, "in")
  )) +
  geom_sf(data = coast, colour = "black", linewidth = 0.1) +
  coord_sf(
    xlim = c(xmin(ext_global), xmax(ext_global)),
    ylim = c(ymin(ext_global), ymax(ext_global)),
    expand = FALSE, clip = "on"
  )

ggsave(
  filename = file.path(project_root, "data/figures/4_01_r_H_TWI_map_mask_with_used_land.png"),
  plot = p_used, width = 14, height = 7, dpi = 600, units = "in"
)

## 3.2. Natural land (fused < 0.05)
p_natural <- plot_cor_twi_vegh(
  input = cor_natural,
  extent = ext_global,
  title = "Pearson’s r (H~TWI) on natural land (fused < 0.05)",
  text_size = 12,
  x_step = 30, y_step = 30
) +
  guides(fill = guide_colorbar(
    title.position = "left",
    barwidth = grid::unit(0.1, "in"),
    barheight = grid::unit(5, "in")
  )) +
  geom_sf(data = coast, colour = "black", linewidth = 0.1) +
  coord_sf(
    xlim = c(xmin(ext_global), xmax(ext_global)),
    ylim = c(ymin(ext_global), ymax(ext_global)),
    expand = FALSE, clip = "on"
  )

ggsave(
  filename = file.path(project_root, "data/figures/4_01_r_H_TWI_map_mask_with_natural_land.png"),
  plot = p_natural, width = 14, height = 7, dpi = 600, units = "in"
)

# ------------------------- 4. Create Distribution Comparison Plot --------------
# Convert rasters to data frames for distribution comparison
df_cor_used <- as.data.frame(cor_used, xy = FALSE, na.rm = TRUE)
df_cor_natural <- as.data.frame(cor_natural, xy = FALSE, na.rm = TRUE)

# Add source labels
df_cor_used$source <- "Used land (fused ≥ 0.05)"
df_cor_natural$source <- "Natural land (fused < 0.05)"

# Combine datasets
df_all <- rbind(df_cor_used, df_cor_natural)
colnames(df_all) <- c("value", "source")

# Define color palette
my_colors <- c(
  "Used land (fused ≥ 0.05)" = "#F8766D",  # red
  "Natural land (fused < 0.05)" = "#00BFC4" # blue
)

# Create distribution comparison plot
text_size <- 12
p_vs <- ggplot(df_all, aes(x = value, fill = source, color = source)) +
  geom_density(alpha = 0.5, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(
    title = NULL,
    x = "Pearson's r(H~TWI)",
    y = "Density",
    fill = "Land type",
    color = "Land type"
  ) +
  scale_fill_manual(values = my_colors, guide = "legend") +
  scale_color_manual(values = my_colors, guide = "legend") +
  theme_bw(base_size = text_size) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = text_size),
    legend.title = element_text(size = text_size),
    axis.title = element_text(size = text_size),
    axis.text = element_text(size = text_size * 0.9),
    plot.title = element_text(size = text_size * 1.1, face = "bold"),
    plot.title.position = "panel"
  )

# ------------------------- 5. Create Scatter Plots by Land-Use Bins ------------
# Stack and convert to data frame for scatter analysis
df <- as.data.frame(c(cor_twi_vegh_r, fused_5km_r, mi_5km_r), xy = TRUE, na.rm = TRUE)
colnames(df) <- c("lon", "lat", "cor", "fused", "mi")

# Create bins for 'fused'
df_binned <- df |>
  filter(!is.na(fused)) |>
  mutate(
    fused_bin = cut(
      fused,
      breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
      include.lowest = TRUE,
      right = FALSE,
      labels = c("0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.2",
                 "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6",
                 "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"),
      ordered_result = TRUE
    )
  )
bins <- levels(df_binned$fused_bin)

# Initialize an empty data frame in long format for regression coefficients
slope_df <- data.frame(
  fused_bin = character(),
  Coefficient = character(),
  Value = numeric(),
  stringsAsFactors = FALSE
)

# Define parameters for layout
x_label <- "Moisture Index (MI)"
y_label <- "Pearson's r(H~TWI)"
n_cols <- 6
n_rows_plots <- ceiling(length(bins) / n_cols)

# Create scatter plots for each bin
plot_list <- list()

for (i in seq_along(bins)) {
  bin_label <- bins[i]
  df_sub <- df_binned |> filter(fused_bin == bin_label)

  # Fit linear model
  lm_fit <- lm(cor ~ mi, data = df_sub)
  slope <- round(coef(lm_fit)[2], 3)
  intercept <- round(coef(lm_fit)[1], 3)

  # Save results in long format
  slope_df <- rbind(
    slope_df,
    data.frame(fused_bin = bin_label, Coefficient = "slope", Value = slope),
    data.frame(fused_bin = bin_label, Coefficient = "intercept", Value = intercept)
  )

  # Create hexbin plot with bin label as title
  p_scatter <- plot_hex_scatter(df_sub, x_var = "mi", y_var = "cor",
                                x_text = NULL, y_text = NULL, title_text = NULL, text_size = 12) +
    labs(x = x_label, y = y_label, title = bin_label) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

  # Adjust axis labels based on position
  row <- ((i - 1) %/% n_cols) + 1
  col <- ((i - 1) %% n_cols) + 1

  if (row < n_rows_plots) {
    p_scatter <- p_scatter + theme(axis.title.x = element_blank())
  }
  if (col > 1) {
    p_scatter <- p_scatter + theme(axis.title.y = element_blank())
  }

  plot_list[[i]] <- p_scatter
}

# Create coefficient plot
slope_df$fused_bin <- factor(slope_df$fused_bin, levels = bins, ordered = TRUE)

p_coe <- ggplot(slope_df, aes(x = fused_bin, y = Value, color = Coefficient, group = Coefficient)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  labs(
    x = "Fused Bin (Land-Use Intensity)",
    y = "Coefficient Value",
    color = "Coefficient"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = text_size),
    legend.title = element_text(size = text_size),
    axis.title = element_text(size = text_size),
    axis.text.x = element_text(angle = 15, hjust = 1),
    axis.text = element_text(size = text_size * 0.9),
    plot.title = element_text(size = text_size * 1.1, face = "bold"),
    plot.title.position = "panel"
  )

# ------------------------- 6. Combine All Visualizations -----------------------
# Manually add an A label to the first subgraph
plot_list[[1]] <- plot_list[[1]] +
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 16),
        plot.tag.position = c(0.05, 0.95))

# Create scatter plots grid
scatter_plots <- wrap_plots(plot_list, ncol = n_cols)

# Combine all components
final_plot <-
  scatter_plots /
  ((p_coe + labs(tag = "B")) | (p_vs + labs(tag = "C"))) +
  plot_layout(heights = c(n_rows_plots, 1)) +
  plot_annotation(
    title = "Land use influence on Pearson's r (H~TWI)",
    subtitle = "Scatter plots across land-use intensity bins (A) with regression coefficients (B) and distribution comparison on natural land and used land (C)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      plot.tag = element_text(size = 16, face = "bold")
    )
  )

# ------------------------- 7. Save All Plots -----------------------------------
# Save combined plot
ggsave(
  filename = here::here("data/figures/4_01_land_use_analysis_r_H_TWI.png"),
  plot = final_plot,
  width = 14,
  height = 10,
  dpi = 600
)

message("✅ All plots saved successfully!")
message("📁 Individual maps saved:")
message("   - 4_01_global_r_H_TWI_on_used_land.png")
message("   - 4_01_global_r_H_TWI_on_natural_land.png")
message("📊 Combined analysis saved:")
message("   - 4_01_combined_analysis_r_H-TWI.png")
message("   - 4_01_r_H_TWI_comparison_on_natural_used_land.png")
