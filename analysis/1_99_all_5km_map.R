# ==============================================================================
# Script: plot_fraction_land_use_maps.R
#
# Purpose:
#   Visualize global fraction land cover layers and annual mean temperature (MAT)
#   using rasters produced from the processing pipeline.
#
# Runtime:
#   ~ 10 mins
# ==============================================================================

# ------------------------- 1. Setup ---------------------------------------------
library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(RColorBrewer)
library(tidyterra)
library(here)

# Load configuration and plotting helper
source(here::here("R/config.R"))
source(here::here("R/plot_var.R"))

message("🌍 Starting visualization maps...")

# ------------------------- 2. Load Data ----------------------------------------
# Load coastlines and land polygons
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
land  <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
land_vect <- vect(land)

# ------------------------- 3. Plot & Save --------------------------------------
plot_map <- function(input, label, output_file) {
  r_layer <- rast(input) |>
    aggregate(c(2,2))

  r_layer <- mask(crop(r_layer, land_vect), land_vect)

  #  r_layer <- mask(crop(r_layer, land_vect), land_vect)
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = r_layer, maxcell = Inf) +
    scale_fill_gradientn(
      colours = rev(brewer.pal(7, "Spectral")),
      na.value = NA) +
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

  ggsave(
    filename = output_file,
    plot     = p,
    width    = 7,
    height   = 3,
    dpi      = 600,
    units    = "in"
  )

  if (file.exists(output_file)) message("✅ Saved map: ", output_file)
}

# Plot maps
plot_map(fused_5km_file, expression(f[used]), here::here("data/figures/1_99_fused_5km_map_0p1d.png"))
plot_map(fbare_5km_file, expression(f[bare]), here::here("data/figures/1_99_fbare_5km_map_0p1d.png"))
plot_map(fwater_5km_file, expression(f[waetr]), here::here("data/figures/1_99_fwater_5km_map_0p1d.png"))
plot_map(fsnow_5km_file, expression(f[snow]), here::here("data/figures/1_99_fsnow_5km_map_0p1d.png"))
# plot_map(mat_1km_file, "Annual mean surface temperature (℃)", here::here("data/figures/1_99_mat_5km_map_0p1d.png"))
# plot_map(map_1km_file, "Annual mean precipitation (mm)", here::here("data/figures/1_99_map_5km_map_0p1d.png"))
# plot_map(srad_1km_file, "Incident solar radiation (kJ m⁻² day⁻¹)", here::here("data/figures/1_99_srad_5km_map_0p1d.png"))
# plot_map(ecoregion_5km_path, "Biomes distribution", here::here("data/figures/1_99_biome_5km_map_0p1d.png"))
