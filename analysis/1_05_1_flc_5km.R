# ==============================================================================
# Purpose:
#   Process global land cover tiles to calculate fractions of:
#     - Vegetation cover (fused)
#     - Bare ground (fbare)
#     - Water (fwater)
#     - Snow/ice (fsnow)
#   at 5km resolution.
#
# Run Time: ~7 mins on UBELIX with 16 core
# ==============================================================================

# -------------------- 1. Setup Environment ------------------------------------
library(terra)
library(dplyr)
library(tidyr)
library(fs)
library(future)
library(furrr)

# Load custom functions
source(here::here("R/config.R"))
source(here::here("R/create_spatial_windows.R"))
source(here::here("R/create_aligned_template.R"))
source(here::here("R/calculate_fraction_land_cover.R"))
source(here::here("R/mosaic_tiles.R"))
source(here::here("R/raster_preprocess_save.R"))

# Set worker numbers for different system
if (hostname == "dash") workers = 4 else workers = 16
message("→ using ", workers, " workers")

# Create output directory
if (!dir.exists(flc_tile_dir)) dir.create(flc_tile_dir, recursive = TRUE)

# # Load Tile Information
tiles_info <- readRDS(valid_tiles_info_path)

t_start <- Sys.time()
message("⏱ Land-cover fraction pipeline started at: ", format(t_start, "%Y-%m-%d %H:%M:%S"))
# -------------------- 3. Parallel Tile Processing -----------------------------
# Single tile processing function
process_tile <- function(tile_info) {

  tile_id <- tile_info$tile_id
  output_file <- file.path(flc_tile_dir, paste0("flc_5km_", tile_id, ".nc"))

  if (file.exists(output_file)) {
    message("⏭️ Skip: ", tile_id)
    return(NULL)
  }

  message("🔹 Processing: ", tile_id)

  ext <- terra::ext(tile_info$xmin, tile_info$xmax, tile_info$ymin, tile_info$ymax)
  lc_r <- terra::rast(cci_landcover_path, lyrs = "lccs_class")
  rc <- terra::crop(lc_r, ext)

  df_win <- create_spatial_windows(rc, value_vars = "lccs_class", dwin = 0.05)
  df_flc <- calculate_fraction_land_cover(df_win, output_file = output_file)

  print(rast(output_file))

  message("✅ Done: ", tile_id)
}

# Parallel processing
run_parallel <- function(tiles_info, n_cores = NULL) {

  if (is.null(n_cores)) n_cores <- availableCores() - 1

  plan(multisession, workers = n_cores)

  dir.create(flc_tile_dir, showWarnings = FALSE, recursive = TRUE)

  future_map(1:nrow(tiles_info), function(i) {
    process_tile(tiles_info[i, ])
  }, .progress = FALSE)

  message("🎉 All tiles completed!")
}

# Usage
run_parallel(tiles_info, n_cores = workers)
# -------------------- 4. Mosaic All Tiles -------------------------------------
message("🗺 Mosaicing all tiles into global raster...")

# Create template raster aligned to 5km grid and mosaic
align_template_5km <- create_aligned_template(twi_450m_mosaic_clean_path, res_out = 0.05)
mosaic_r <- mosaic_tiles(
  input_dir   = flc_tile_dir,
  output_file = NULL,
  pattern = "*.nc",
  target_grid = align_template_5km,
  if_resample = FALSE,
  if_crop = TRUE
)

print(rast(mosaic_r))
print(summary(mosaic_r))

# -------------------- 5. Save to Single Layer Raster -------------------------
message("🔧 Saving single-layer rasters...")

# Set output files and variable names
output_files <- c(fused_5km_file, fbare_5km_file, fwater_5km_file, fsnow_5km_file)
varnames     <- c("fused", "fbare", "fwater", "fsnow")

# save rasters to single-layer rasters
raster_preprocess_save(
  input            = mosaic_r,
  output           = output_files,
  varname          = varnames,
  if_aggregate     = FALSE,
  if_resample      = FALSE,
  if_return_raster = FALSE
)

elapsed_total <- difftime(Sys.time(), t_start, units = "mins")
message(sprintf("🎉 Pipeline completed successfully [%.1f mins]", elapsed_total))

print(rast(fused_5km_file))
# # -------------------- 6. Cleanup Intermediate Tiles ---------------------------
# tiles_path <- fs::dir_ls(path = flc_tile_dir, glob = "*.nc")
# if (length(tiles_path) > 0) {
#   message("🧹 Cleaning intermediate tile files...")
#   file.remove(tiles_path)
#   message("✅ Intermediate files removed.")
# }

# -------------------- 7. Optional Check ---------------------------
# r <- rast(fused_5km_file)
# print(r)
# summary(r)
# plot(r, main = "Used land fraction (5km)")
#
# r <- rast("/storage/scratch/giub_geco/tting/data/global_flc_5km/30_30_deg/flc_60S_90W.nc")
# print(r)
# summary(r)
# plot(r, main = "Bare land fraction (5km)")
#
# r <- rast(fwater_5km_file)
# print(r)
# summary(r)
# plot(r, main = "Water body fraction (5km)")
#
# r <- rast(fsnow_5km_file)
# print(r)
# summary(r)
# plot(r, main = "Permanent snow and ice fraction (5km)")
#
# # single tile check
# args <- tiles_info[11,]
#
# tile_id <- args$tile_id
# output_file <- file.path(flc_tile_dir, paste0("flc_5km_", tile_id, ".nc"))
#
# r_out <- rast(output_file)
# print(r_out)
