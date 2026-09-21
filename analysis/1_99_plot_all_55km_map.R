# ==============================================================================
# Script: plot_environmental_maps_55km.R
#
# Purpose:
#   Visualize multiple environmental variables at 55km resolution, including:
#       - Moisture Index (MI)
#       - Fraction land use (fused, fbare, fwater, fsnow)
#       - Elevation metrics (dem_sd, dem_rg)
#       - Protected area fraction (fpa)
#       - Annual mean air temperature (mat)
#
# Output:
#   - High-resolution PNG maps saved in data/figures/
# ==============================================================================

# ------------------------- 1. Setup Environment --------------------------------
library(terra)
library(tidyterra)
library(purrr)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(RColorBrewer)
library(khroma)
library(here)

# Load configuration and plotting helper
source(here::here("R/config.R"))
source(here::here("R/plot_var.R"))

message("🌍 Starting visualization of 55-km environmental maps...")

# ------------------------- 2. Load Data ----------------------------------------
# Load coastlines and land polygons
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")
land  <- rnaturalearth::ne_countries(scale = 110, returnclass = "sf")
land_vect <- vect(land)

# ------------------------- 3. Plot & Save --------------------------------------
plot_map <- function(input, title_text, output_file) {
  r_layer <- rast(input)
  #  r_layer <- mask(crop(r_layer, land_vect), land_vect)
  p <- plot_var(
    input      = r_layer,
    title_text = title_text,
    extent     = ext_global,
    text_size  = 12,
    x_step     = 30,
    y_step     = 30
  ) +
    guides(fill = guide_colorbar(
      title.position = "left",
      barwidth      = grid::unit(0.1, "in"),
      barheight     = grid::unit(5, "in")
    )) +
    geom_sf(data = coast, colour = "black", linewidth = 0.1) +
    coord_sf(
      xlim = c(terra::xmin(ext_global), terra::xmax(ext_global)),
      ylim = c(terra::ymin(ext_global), terra::ymax(ext_global)),
      expand = FALSE,
      clip = "on"
    )

  ggsave(
    filename = output_file,
    plot     = p,
    width    = 14,
    height   = 7,
    dpi      = 300,
    units    = "in"
  )

  if (file.exists(output_file)) message("✅ Saved map: ", output_file)
}

# Plot each environmental variable
plot_map(mi_55km_file,     "Moisture Index", here::here("data/figures/1_99_region_grid_mi_55km_map.png"))
plot_map(fused_55km_file,  "Land use fraction", here::here("data/figures/1_99_region_grid_fused_55km_map.png"))
plot_map(fbare_55km_file,  "Bare land fraction", here::here("data/figures/1_99_region_grid_fbare_55km_map.png"))
plot_map(fwater_55km_file, "Water body fraction", here::here("data/figures/1_99_region_grid_fwater_55km_map.png"))
plot_map(fsnow_55km_file, "Permanent snow & ice fraction", here::here("data/figures/1_99_region_grid_fsnow_55km_map.png"))
plot_map(dem_sd_55km_path,"Elevation SD", here::here("data/figures/1_99_region_grid_dem_sd_55km_map.png"))
plot_map(fpa_55km_path,   "Protected area fraction", here::here("data/figures/1_99_region_grid_fpa_55km_map.png"))
plot_map(dem_rg_95p_05p_55km_path, "Elevation range (95p-05p)", here::here("data/figures/1_99_region_grid_dem_rg_95p_05p_55km_map.png"))
plot_map(fveg_55km_path,   "Vegetated area fraction", here::here("data/figures/1_99_region_grid_fveg_55km_map.png"))
plot_map(fveg_real_55km_path,  "Real vegetated area fraction", here::here("data/figures/1_99_region_grid_fveg_real_55km_map.png"))

# ------------------------- 4. Completion Message --------------------------------
message("🎉 All 55-km environmental maps have been created and saved successfully!")

