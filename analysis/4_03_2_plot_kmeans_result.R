# ========================================================
#  Script: K-means 8-cluster analysis and biome composition
#  Purpose:
#    - Load raster datasets
#    - Summarize K-means (K=8) clusters
#    - Visualize spatial patterns and statistical distributions
#    - Quantify absolute biome composition per cluster
# ========================================================

# ========================================================
# 1. Load required packages
# ========================================================

library(terra)        # raster processing
library(dplyr)        # data manipulation
library(ggplot2)      # plotting
library(patchwork)    # multi-panel layouts
library(tidyterra)    # terra + ggplot2 bridge
library(rnaturalearth)# coastline data
library(sf)           # vector data handling

source(here::here("R/config.R"))
source(here::here("R/plot_boxplot.R"))
source(here::here("R/plot_kmeans_map.R"))

# ========================================================
# 3. Load raster datasets
# ========================================================

# Environmental variables and cluster map (5 km resolution)
mi_5km_r           <- terra::rast(mi_5km_file) * 0.0001
fused_5km_r        <- terra::rast(fused_5km_file)
cor_twi_vegh_5km_r <- terra::rast(cor_twi_vegh_mosaic_file)

# K-means (K = 8) result and biome classification
kmeans_8c_r  <- terra::rast(kmeans_map_8c_path)
ecoregion_r  <- terra::rast(ecoregion_5km_path)

# Quick sanity checks
mi_5km_r
fused_5km_r
cor_twi_vegh_5km_r
kmeans_8c_r
ecoregion_r

# ========================================================
# 4. Raster stacking and dataframe conversion
# ========================================================

# Stack all rasters into a single object
stacked <- c(
  cor_twi_vegh_5km_r,
  fused_5km_r,
  mi_5km_r,
  kmeans_8c_r,
  ecoregion_r
)

# Convert raster stack to point-wise dataframe
df <- as.data.frame(stacked, xy = TRUE, na.rm = TRUE)

# Explicit column naming for downstream clarity
colnames(df) <- c(
  "lon", "lat",
  "cor", "fused", "mi",
  "cluster8c",
  "BIOME_NUM"
)

# Free memory (large rasters no longer needed)
rm(
  mi_5km_r,
  fused_5km_r,
  cor_twi_vegh_5km_r,
  ecoregion_r,
  stacked
)
gc()

# ========================================================
# 5. Cluster-level statistics
# ========================================================

# Summarize distribution of key variables per cluster
df_summary <- df |>
  group_by(cluster8c) |>
  summarise(
    count      = n(),
    percentage = n() / nrow(df) * 100,

    Q1_mi      = quantile(mi, 0.25),
    median_mi  = median(mi),
    Q3_mi      = quantile(mi, 0.75),

    Q1_cor     = quantile(cor, 0.25),
    median_cor = median(cor),
    Q3_cor     = quantile(cor, 0.75),

    Q1_fused   = quantile(fused, 0.25),
    median_fused = median(fused),
    Q3_fused   = quantile(fused, 0.75),

    .groups = "drop"
  ) |>
  mutate(percentage = round(percentage, 2)) |>
  arrange(median_mi)   # order clusters along MI gradient

message("df_summary:")
print(df_summary$median_mi)
print(df_summary$cluster8c)
print(df_summary$percentage)

# ========================================================
# 6. Assign descriptive cluster labels and colors
# ========================================================

# Cluster IDs ordered by median MI
cluster_values <- df_summary$cluster8c

# Human-readable labels (interpretation-driven)
cluster_labels <- c(
  "Arid downslope \nsuppression",
  "Arid downslope \nsupport",
  "Upslope \nland use",
  "Transition",
  "Downslope \nland use",
  "Moist downslope \nsupression",
  "Moist downslope \nsupport",
  "Super-moist \ndownslope suppression"
)

# Map labels to original cluster IDs
names(cluster_labels) <- cluster_values

# Fixed color palette for all figures
fill_colors <- setNames(
  c(
    "#E78AC3", "#FC8D62", "#FFD92F", "#E5C494",
    "#B3B3B3", "#66C2A5", "#8DA0CB", "#A6D854"
  ),
  cluster_labels
)

# Apply ordered factor with descriptive labels
df$cluster8c <- factor(
  df$cluster8c,
  levels = cluster_values,
  labels = cluster_labels
)

# Apply the same factor ordering to summary table
df_summary <- df_summary |>
  mutate(
    cluster8c = factor(
      cluster_labels[as.character(cluster8c)],
      levels = cluster_labels
    )
  )

# ========================================================
# 7. Global K-means cluster map
# ========================================================

text_size <- 7

# Load global coastline for map context
coast <- rnaturalearth::ne_coastline(
  scale = 110,
  returnclass = "sf"
)

# Base global cluster map
p_8c <- plot_kmeans_map(
  kmeans_8c_r,
  text_size  = text_size,
  extent     = ext_global,
  title_text = ""
) +
  geom_sf(data = coast, colour = "black", linewidth = 0.1) +
  coord_sf(
    xlim   = c(xmin(ext_global), xmax(ext_global)),
    ylim   = c(ymin(ext_global), ymax(ext_global)),
    expand = FALSE,
    clip   = "on"
  ) +
  theme(
    legend.position  = "none",
    axis.title       = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_blank()
  ) +
  labs(tag = "a)") +
  theme(
    plot.tag = element_text(size = text_size, face = "bold"),
    plot.tag.position = c(0.01, 1),
    plot.title = element_blank()
  )

# ========================================================
# 8. Cluster-wise variable distributions
# ========================================================

# MI distribution

pbox_mi <- plot_boxplot(
  df, "cluster8c", "mi",
  expression(MI),
  text_size = text_size,
  show_legend = FALSE,
  ylim = c(0,3.5)
) +
  labs(tag = "b)") +
  theme(
    plot.tag = element_text(size = text_size, face = "bold"),
    plot.tag.position = c(0.03, 1)
  )

# Fused index distribution
pbox_fused <- plot_boxplot(
  df, "cluster8c", "fused",
  bquote(f[.("used")]),
  text_size = text_size,
  show_legend = FALSE
) + labs(tag = "c)") +
  theme(
    plot.tag = element_text(size = text_size, face = "bold"),
    plot.tag.position =  c(0.03, 1)
  )

# Correlation distribution
pbox_cor <- plot_boxplot(
  df, "cluster8c", "cor",
  bquote(r[.("H, TWI")]),
  text_size = text_size,
  show_legend = FALSE,
  ylim = c(-0.75,0.75)
) + labs(tag = "d)") +
  theme(
    plot.tag = element_text(size = text_size, face = "bold"),
    plot.tag.position =  c(0.03, 1)
  )

# ========================================================
# 9. Cluster area percentage
# ========================================================

# Percentage of global pixels per cluster
p_bar <- ggplot(
  df_summary,
  aes(x = cluster8c, y = percentage, fill = cluster8c)
) +
  geom_col(width = 0.9, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = sprintf("%.1f", percentage)),
    vjust = -0.5,
    size = 1.5
  ) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_x_discrete(drop = TRUE, expand = c(0.1, 0.1)) +
  theme_bw(base_size = text_size) +
  theme(
    legend.position = "none",
    axis.text.x     = element_blank()
  ) +
  labs(
    x = "",
    y   = "Frequency (%)",
    tag = "e)"
  ) +
  theme(
    plot.tag = element_text(size = text_size, face = "bold"),
    plot.tag.position = c(0.03, 1)
  )


# ========================================================
# 10. Combine map and statistics into one figure
# ========================================================

# Extract legend from map
p_8c_legend <- p_8c + theme(legend.position = "bottom")

# Final multi-panel layout
p_8c_stat <- wrap_plots(
  p_8c,
  wrap_plots(pbox_mi, pbox_fused, pbox_cor, p_bar, nrow = 1),
  wrap_elements(cowplot::get_legend(p_8c_legend)),
  nrow = 3,
  heights = c(3.2, 1.3, 0.5)
) &
  theme(
    plot.margin = margin(1, 0, 0, 1, "pt"),  # 移除所有子图边距
    plot.background = element_blank()
  )

# Save high-resolution figure
ggsave(
  filename = here::here("data/figures/4_03_kmeans_gl_map_8c_stat.png"),
  plot     = p_8c_stat,
  width    = 7,
  height   = 5,
  dpi      = 600,
  units    = "in"
)


# ========================================================
# 11. Individual cluster maps (1–8)
# ========================================================

plot_list <- list()

for (i in seq_along(cluster_labels)) {

  cluster <- gsub("\n", "", cluster_labels[i])

  p <- plot_kmeans_map(
    kmeans_8c_r,
    text_size = text_size,
    extent = ext_global,
    title_text = paste0(cluster),
    highlight_cluster = cluster_values[i]
  ) +
    geom_sf(data = coast, colour = "black", linewidth = 0.1) +
    coord_sf(
      xlim   = c(xmin(ext_global), xmax(ext_global)),
      ylim   = c(ymin(ext_global), ymax(ext_global)),
      expand = FALSE
    ) +
    theme(
      legend.position = "none",
      axis.title      = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(tag = paste0(letters[i], ")"))

  plot_list[[i]] <- p
}

p_legend <- plot_list[[1]] + theme(legend.position = "bottom")

p_1_to_8c <- wrap_plots(plot_list, ncol = 2) /
  wrap_elements(cowplot::get_legend(p_legend)) +
  plot_layout(heights = c(5, 0.3))

ggsave(
  filename = here::here("data/figures/4_03_kmeans_gl_map_1_to_8.png"),
  plot     = p_1_to_8c,
  width    = 7,
  height   = 7.5,
  dpi      = 600,
  units    = "in"
)

# ========================================================
# 12. Biome composition (absolute counts)
# ========================================================

ecoregion <- vect(ecoregion_path)

biomes_info <- ecoregion |>
  as.data.frame() |>
  select(BIOME_NUM, BIOME_NAME, COLOR_BIO) |>
  distinct() |>
  arrange(BIOME_NUM)

df_biome_summary_counts <- df |>
  group_by(cluster8c, BIOME_NUM) |>
  summarise(count = n(), .groups = "drop") |>
  left_join(
    biomes_info |> select(BIOME_NUM, BIOME_NAME, COLOR_BIO),
    by = "BIOME_NUM"
  ) |>
  mutate(cluster8c = factor(cluster8c, levels = levels(df$cluster8c))) |>
  group_by(cluster8c) |>
  arrange(cluster8c, desc(count)) |>
  mutate(BIOME_NAME = factor(BIOME_NAME, levels = unique(BIOME_NAME))) |>
  ungroup()


# ========================================================
# 13. Plot absolute biome composition
# ========================================================

p_8c_biome_counts <- ggplot(
  df_biome_summary_counts,
  aes(x = cluster8c, y = count, fill = BIOME_NAME)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = setNames(biomes_info$COLOR_BIO, biomes_info$BIOME_NAME),
    name = "Biome"
  ) +
  scale_y_continuous(
    name   = "Number of Observations (×10⁴)",
    labels = function(x) x / 10000
  ) +
  labs(
    x     = "Cluster",
    title = "Absolute Biome Composition of Each Cluster"
  ) +
  theme_bw(base_size = text_size) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = text_size * 1.2, face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 2))

ggsave(
  filename = here::here("data/figures/4_03_kmeans_8c_biome_counts.png"),
  plot     = p_8c_biome_counts,
  width    = 7,
  height   = 7,
  dpi      = 300
)
