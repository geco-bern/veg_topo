# UBELIX: 20 cores, 400G, done [7-10 min]

# ------Load required libraries-------------------------------------------------

library(terra)
library(dplyr)
library(furrr)

source(here::here("R/config.R"))
source(here::here("R/mosaic_tiles.R"))
source(here::here("R/extent_to_tile_ids.R"))
source(here::here("R/create_spatial_windows.R"))
source(here::here("R/calculate_correlation_bywin.R"))
source(here::here("R/raster_preprocess_save.R"))

# Set worker numbers for different system
hostname <- trimws(tolower(system("hostname", intern = TRUE)))
if (hostname == "dash") {
  workers = 2
  message("→ using ", workers, " workers")

} else {
  workers = 17
  message("→ using ", workers, " workers")
}

# ------ File Configuration ---------------------------------------------

if (!dir.exists(cor_twi_vegh_tiles_dir)) {
  dir.create(cor_twi_vegh_tiles_dir, recursive = TRUE)
  message("Directory created: ", cor_twi_vegh_tiles_dir)
}

# ------------ function to process single tile----------------------------------
process_r_H_TWI_5km <- function(tile_row, output_dir = cor_twi_vegh_tiles_dir) {

  tryCatch({

    # --- Tile info ---
    tile_id <- tile_row$tile_id
    tile_extent <- terra::ext(tile_row$xmin, tile_row$xmax, tile_row$ymin, tile_row$ymax)
    tile_xmid <- (tile_row$xmin + tile_row$xmax)/2
    tile_ymid <- (tile_row$ymin + tile_row$ymax)/2

    tictoc::tic(paste0("Processing tile: ", tile_id))
    t0 <- Sys.time()

    # --- TWI Raster ---
    twi_rc <- terra::rast(twi_450m_mosaic_clean_path) |> terra::crop(tile_extent)
    twi_tmp_path <- file.path(tempdir(), paste0("tile_", tile_id, "_twi_450m.nc"))
    terra::writeCDF(twi_rc, twi_tmp_path, varnames = "twi", overwrite = TRUE)
    rm(twi_rc); gc()
    twi_rc <- terra::rast(twi_tmp_path)
    names(twi_rc) <- "twi"
    message("Saved temporary TWI raster: ", twi_tmp_path)

    # --- Vegetation Height Raster ---
    vegh_rc <- terra::rast(vegh_450m_mosaic_path) |> terra::crop(tile_extent)
    vegh_tmp_path <- file.path(tempdir(), paste0("tile_", tile_id, "_vegh_450m.nc"))
    terra::writeCDF(vegh_rc, vegh_tmp_path, varnames = "vegh",overwrite = TRUE)
    rm(vegh_rc); gc()
    vegh_rc <- terra::rast(vegh_tmp_path)
    names(vegh_rc) <- "vegh"
    message("Saved temporary vegetation height raster: ", vegh_tmp_path)

    # --- Stack and correlation ---
    stacked <- c(twi_rc, vegh_rc)
    df_win <- create_spatial_windows(stacked, value_vars = c("twi", "vegh"), dwin = 0.05)
    df_cor <- calculate_correlation_bywin(df_win, x = "twi", y = "vegh")

    # --- Save correlation ---
    cor_r <- terra::rast(df_cor[, c("lon_mid", "lat_mid", "correlation")], type="xyz", crs="EPSG:4326")
    names(cor_r) <- "r_H_TWI"
    cor_nc_path <- file.path(output_dir, paste0("tile_", tile_id, "_r_H_TWI_5km_map.nc"))
    terra::writeCDF(cor_r, cor_nc_path, overwrite = TRUE)
    message("Saved: ", cor_nc_path)

    # --- Save p-value ---
    pval_r <- terra::rast(df_cor[, c("lon_mid", "lat_mid", "cor_pval")], type="xyz", crs="EPSG:4326")
    names(pval_r) <- "pval_r_H_TWI"
    pval_nc_path <- file.path(output_dir, paste0("tile_", tile_id, "_r_H_TWI_5km_pval.nc"))
    terra::writeCDF(pval_r, pval_nc_path, overwrite = TRUE)
    message("Saved: ", pval_nc_path)

    # --- Print proccess time ---
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("Region %s completed [%.1f mins]", tile_id, elapsed_mins))
    tictoc::toc()

    return(TRUE)

  }, error = function(e) {
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("❌ Tile %s failed after %.1f mins: %s", tile_row$tile_id, elapsed_mins, e$message))
    return(FALSE)
  })
}

# # ----------------- Parallel execution for all tiles-----------------
# # Set up cluster plan
# plan(cluster, workers = workers)
#
# # load tiles info
# tiles_info <- readRDS(valid_tiles_info_path)
#
# # Run in parallel
# tictoc::tic("🚀 Parallel processing of tiles")
# results <- future_map(seq_len(nrow(tiles_info)),
#                       function(i) process_r_H_TWI_5km(tiles_info[i, ]),
#                       .progress = FALSE,
#                       .options = furrr_options(seed=TRUE))
# plan(sequential)
# tictoc::toc()
#
# # Summarize results
# success_count <- sum(unlist(results))
# fail_count <- length(results) - success_count
# message(sprintf("✅ Completed: %d succeeded, ❌ %d failed.", success_count, fail_count))

# -------- Combination ---------------------------------------------------------

# mosacing the r(H~TWI) map
mosaic_tiles(
  input_dir   = cor_twi_vegh_tiles_dir,
  output_file = cor_twi_vegh_mosaic_file,
  pattern = "*_map.nc",
  varname = "r_H_TWI")

# mosacing the pval ofr(H~TWI) map
mosaic_tiles(
  input_dir   = cor_twi_vegh_tiles_dir,
  output_file = pval_cor_twi_vegh_mosaic_file,
  pattern = "*_pval.nc",
  varname = "p_H_TWI")

# # ---------- Delete intermediate data ------------------------------------------
# # List all files in the directory cor_twi_vegh_tiles_dir that match "*.nc"
# # If there are any files found, delete them
#
# cor_5km_tiles_path <- fs::dir_ls(path = cor_twi_vegh_tiles_dir, glob = "*.nc")
# if (length(cor_5km_tiles_path) > 0) file.remove(cor_5km_tiles_path)

# ------ Single tile test ---------------------------------------------

# test for 1 tile
# tiles_info <- readRDS(valid_tiles_info_path)
# process_r_H_TWI_5km(tiles_info[35,])

# ------ Smaller region  test ---------------------------------------------
# regA_info <- data.frame(
#   tile_id = c("Aletsch_glacier"),
#   ymin = c(46.8),
#   ymax = c(47),
#   xmin = c(7.8),
#   xmax = c(8),
#   sample_id = c(1)
# )
#
# # center location
# regA_info$xmid <- (regA_info$xmax + regA_info$xmin) / 2
# regA_info$ymid <- (regA_info$ymax + regA_info$ymin) / 2
#
# tile_row <- regA_info[1, ]
#
# output_dir = cor_twi_vegh_tiles_dir
# text_size = 12
# fig_width = 14
# fig_height = 14



