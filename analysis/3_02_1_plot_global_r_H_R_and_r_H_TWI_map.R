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
# 1. Data preparation function
# =============================================================================

prepare_correlation_data <- function(raster_obj, var_name) {
  # Convert to data frame and clean
  df <- as.data.frame(raster_obj, xy = FALSE) |>
    drop_na()
  colnames(df)[1] <- "correlation"

  # Create grouping variable
  df <- df |>
    mutate(group = case_when(
      correlation > 0 ~ "Positive",
      correlation < 0 ~ "Negative"
    ))

  # Calculate percentages
  percent_positive <- sum(df$correlation > 0) / nrow(df) * 100
  percent_negative <- sum(df$correlation < 0) / nrow(df) * 100

  return(list(
    data = df,
    percent_positive = percent_positive,
    percent_negative = percent_negative
  ))
}

create_correlation_histogram <- function(data_list, var_name, colors, y_max = NULL) {
  # Extract colors
  positive_color <- colors[50]
  negative_color <- colors[200]

  # Get max count for y-axis if not provided
  if (is.null(y_max)) {
    bin_counts <- table(cut(data_list$data$correlation, breaks = seq(-1, 1, by = 0.05)))
    y_max <- max(bin_counts) * 1.15
  }

  p <- ggplot(data_list$data, aes(x = correlation, fill = group)) +
    geom_histogram(binwidth = 0.05, color = "black", linewidth = 0.2,
                   alpha = 0.5, boundary = 0) +
    scale_fill_manual(
      values = c("Positive" = positive_color, "Negative" = negative_color),
      name = "Correlation Type"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
    xlim(-1, 1) +
    labs(
      title = NULL,
      x = var_name,
      y = expression("Count (×10"^4*")")
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      labels = function(x) round(x / 10000, 1)
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
      axis.title.x = element_text(hjust = 0.5, margin = margin(t = 1), size = 7),
      axis.title.y = element_text(hjust = 0.5, angle = 90, margin = margin(r = 1), size = 7),
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.ticks.x = element_line(linewidth = 0.3),
      axis.ticks.y = element_line(linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.2),
      panel.background = element_blank()
    ) +
    annotate("text", x = 0.5, y = y_max * 0.85,
             label = sprintf("%.1f%%", data_list$percent_positive),
             color = positive_color, size = 2.5,
             fontface = "bold",
             vjust = -0.5) +
    annotate("text", x = -0.5, y = y_max * 0.85,
             label = sprintf("%.1f%%", data_list$percent_negative),
             color = negative_color, size = 2.5,
             fontface = "bold",
             vjust = -0.5)

  return(p)
}

create_correlation_map <- function(raster_obj, palette_name, label) {
  # Aggregate for coarser resolution
  raster_agg <- aggregate(raster_obj, fact = c(2, 2))

  # Load coastline
  coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = raster_agg, maxcell = Inf) +
    scico::scale_fill_scico(
      palette = palette_name,
      direction = -1,
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.5),
      midpoint = 0,
      na.value = NA
    ) +
    ggplot2::labs(
      title = NULL,
      fill = label,
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(from = -180, to = 180, by = 30),
      limits = c(-180, 180),
      expand = expansion(mult = 0.0001)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(from = -60, to = 90, by = 30),
      limits = c(-60, 90),
      expand = expansion(mult = 0.0001)
    ) +
    ggplot2::theme_bw(base_size = 7) +
    guides(fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(0.1, "in"),
      barheight = grid::unit(2.6, "in")
    )) +
    geom_sf(data = coast,
            colour = 'black',
            linewidth = 0.1) +
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

  return(p)
}

# =============================================================================
# 2. Load and prepare data for both variables
# =============================================================================

# 2.1 Load coastline data
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# 2.2 Read correlation rasters
r_H_Rin <- rast(r_H_R_5km_path)
r_H_TWI <- rast(cor_twi_vegh_mosaic_file)

# 2.3 Prepare data for Rin
data_rin <- prepare_correlation_data(r_H_Rin, "r(H,Rᵢₙ)")
colors_rin <- scico(256, palette = "vik")

# 2.4 Prepare data for TWI
data_twi <- prepare_correlation_data(r_H_TWI, "r(H,TWI)")
colors_twi <- scico(256, palette = "bam")

# 2.5 Calculate common y-axis limits for histograms
bin_counts_rin <- table(cut(data_rin$data$correlation, breaks = seq(-1, 1, by = 0.05)))
bin_counts_twi <- table(cut(data_twi$data$correlation, breaks = seq(-1, 1, by = 0.05)))
y_max_common <- max(max(bin_counts_rin), max(bin_counts_twi)) * 1.15

# =============================================================================
# 3. Create plots for Rin
# =============================================================================

# 3.1 Histogram
p_hist_rin <- create_correlation_histogram(
  data_rin,
  "r(H,Rᵢₙ)",
  colors_rin,
  y_max = y_max_common
)

# 3.2 Map
p_map_rin <- create_correlation_map(r_H_Rin, "vik", expression(r(H, R["in"]))) +
  labs(tag = "b)") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_blank()
  )

# 3.3 Combine Rin: map + inset histogram
p_comb_rin <- p_map_rin +
  inset_element(
    p_hist_rin +
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

# =============================================================================
# 4. Create plots for TWI
# =============================================================================

# 4.1 Histogram
p_hist_twi <- create_correlation_histogram(
  data_twi,
  "r(H,TWI)",
  colors_twi,
  y_max = y_max_common
)

# 4.2 Map
p_map_twi <- create_correlation_map(r_H_TWI, "bam", "r(H,TWI)") +
  labs(tag = "a)") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_blank()
  )

# 4.3 Combine TWI: map + inset histogram
p_comb_twi <- p_map_twi +
  inset_element(
    p_hist_twi +
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

# =============================================================================
# 5. Combine both panels into final figure
# =============================================================================

# 5.1 Combine using patchwork
p_final <- (p_comb_twi / p_comb_rin) +
  plot_annotation(
    title = NULL,
    theme = theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_blank()
    )
  ) &
  theme(
    plot.background = element_blank()
  )

# 5.2 Save combined figure
ggsave(
  filename = file.path(project_root, "data/figures/3_01_combined_r_H_Rin_TWI_5km_map.png"),
  plot = p_final,
  width = 7,
  height = 6,
  dpi = 600,
  units = "in"
)


cat("\n✅ Done! Saved to: data/figures/3_01_combined_r_H_Rin_TWI_5km_map.png\n")


# =============================================================================
# 5. Plot p-value map (p_pval)
# =============================================================================

# 5.1 Aggregate p-value raster
p_H_Rin <- rast(pval_r_H_R_5km_path) |>
  aggregate(fact = c(2, 2))

# 5.2 Create p-value map using custom function
p_pval_Rin <- plot_cor_pval(
  input = p_H_Rin,
  extent = ext_global,
  title_text = NULL,
  fill_text = "p(H,Rᵢₙ)",
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
  ) +
  labs(tag = "b)") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(5, 0, 0, 0),
    plot.background = element_blank()
  )



# 5.1 Aggregate p-value raster
p_H_TWI <- rast(pval_cor_twi_vegh_mosaic_file) |>
  aggregate(fact = c(2, 2))

# 5.2 Create p-value map using custom function
p_pval_TWI <- plot_cor_pval(
  input = p_H_TWI,
  extent = ext_global,
  title_text = NULL,
  fill_text = "p(H,TWI)",
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
  ) +
  labs(tag = "a)") +
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(0, 0, 5, 0),
    plot.background = element_blank()
  )


# 5.1 Combine using patchwork
p_final_pval <- (p_pval_TWI / p_pval_Rin) +
  plot_annotation(
    title = NULL,
    theme = theme(
      plot.background = element_blank()
    )
  ) &
  theme(
    plot.background = element_blank()
  )

# 5.2 Save combined figure
ggsave(
  filename = file.path(project_root, "data/figures/3_01_combined_p_H_Rin_TWI_5km_map.png"),
  plot = p_final_pval,
  width = 7,
  height = 6,
  dpi = 600,
  units = "in"
)
