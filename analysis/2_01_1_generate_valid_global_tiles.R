# ～1.5 min

# ------Load required libraries-------------------------------------------------
library(terra)
library(sf)
library(parallel)

# ------Load configuration and helper functions---------------------------------
source(here::here("R/config.R"))
source(here::here("R/generate_tile_grid.R"))
source(here::here("R/filter_land_tiles_parallel.R"))

# Set worker numbers for different system
hostname <- trimws(tolower(system("hostname", intern = TRUE)))
if (hostname == "dash") {
  workers = 8
  message("→ using ", workers, " workers")

} else {
  workers = 50
  message("→ using ", workers, " workers")
}
# --------- Generate and valid the land tile with TWI data ---------------------

tictoc::tic()
tile_grid <- generate_tile_grid()
valid_tiles <- filter_land_tiles_parallel(tile_grid, twi_450m_mosaic_clean_path, workers)
saveRDS(valid_tiles, valid_tiles_info_path)
tictoc::toc()

# -------- Create GeoPackage for each valid land tile---------------------------

# Create list of rectangle polygons with tile_id as attribute
rects_list <- lapply(1:nrow(valid_tiles), function(i) {
  # Define polygon geometry using corner coordinates
  geom <- st_sfc(st_polygon(list(rbind(
    c(valid_tiles$xmin[i], valid_tiles$ymin[i]),
    c(valid_tiles$xmin[i], valid_tiles$ymax[i]),
    c(valid_tiles$xmax[i], valid_tiles$ymax[i]),
    c(valid_tiles$xmax[i], valid_tiles$ymin[i]),
    c(valid_tiles$xmin[i], valid_tiles$ymin[i])
  ))), crs = 4326)  # Set CRS to WGS84 (EPSG:4326)

  # Create a single-feature sf object with tile_id as a field
  st_sf(name = valid_tiles$tile_id[i], geometry = geom)
})

# Combine and save the resulting sf object to disk
rects_sf <- do.call(rbind, rects_list)
st_write(rects_sf, valid_tiles_vect_path, layer = "valid_tile", delete_layer = TRUE)

# -------- check valid land tile---------------------------
# tiles_check <- st_read(valid_tiles_vect_path, layer = "valid_tile")
# plot(st_geometry(tiles_check))
