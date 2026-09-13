# ～ 12h on UBELIX
# ❌ Error in validation for "35WNT": [mosaic] resolution does not match

library(terra)
library(dplyr)
library(sf)
library(tidyr)
library(furrr)
library(future)
library(ggplot2)
library(patchwork)

source(here::here("R/config.R"))
source(here::here("R/create_aligned_template.R")) # with source(here::here("R/get_lonlat_extent.R"))
source(here::here("R/create_spatial_windows.R"))
source(here::here("R/calculate_correlation_bywin.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/df_to_raster.R")) # with source(here::here("R/save_raster.R"))
source(here::here("R/aggregate_topography.R"))
source(here::here("R/extent_to_tile_ids.R"))
source(here::here("R/calc_meteoland_sw_in.R"))
source(here::here("R/calc_tri.R"))


# Set worker configuration based on system
if (hostname == "dash") {
  chunk_size <- 1000
  workers <- 16
} else {
  chunk_size <- 10000
  workers <- 100
}
message("→ Using ", workers, " workers and chunk_size = ", chunk_size, " for parallel processing.\n")

# All ALS_MAX data
file_bases <- c("30SWH","30TWN","30UVD","31UFT","32TMS",
                "32TMT", "34WFS", "35VMF","35VNL","35WMR","35WNT",
                "06VXR", "08WNA", "12VUN", "16PHS", "32MPC", "32MQE",
                "06WVT", "10TER", "13UEA", "16SGE", "32MPE", "32NNF")

# Function to process for each Lidar tile
validate_vegetation_heights <- function(file_base) {

  # ============================================================================
  # 1. INITIALIZATION
  # ============================================================================
  message("⭐️⭐️⭐️ Processing:", file_base, " ⭐️⭐️⭐️")

  # Create output directory
  if (!dir.exists(h_validation_dir)) dir.create(h_validation_dir, recursive = TRUE)

  # ----- Check if input file exists -----
  # LiDAR input path and alignment grid
  lidar_path <- file.path(lidar_asl_dir, paste0(file_base, ".tif"))
  lidar_path <- if (file.exists(lidar_path)) lidar_path else file.path(lidar_lvis_dir, paste0(file_base, ".tif"))
  if (!file.exists(lidar_path)) return(message("❌ LiDAR file not found"))

  # ----- Region info -----
  # Create alignment template (30m resolution)
  grid_30m <- create_aligned_template( lidar_path , res_out = 0.00025, trim_input = TRUE)
  reg_extent <- ext(grid_30m)
  message("Extent:",reg_extent[1], ", ", reg_extent[2],", ",reg_extent[3],", ",reg_extent[4], " (xmin,xmax,ymin,ymax)")

  # ----- Output file paths -----
  hlidar_file <- file.path(h_validation_dir, paste0(file_base, "_hlidar.tif"))
  hlang_file <- file.path(h_validation_dir, paste0(file_base, "_hlang.tif"))
  twi_file <- file.path(h_validation_dir, paste0(file_base, "_twi.tif"))
  dem_file <- file.path(h_validation_dir, paste0(file_base, "_dem.tif"))
  slope_file <- file.path(h_validation_dir, paste0(file_base, "_slope.tif"))
  aspect_file <- file.path(h_validation_dir, paste0(file_base, "_aspect.tif"))
  tri_file <- file.path(h_validation_dir, paste0(file_base, "_tri.tif"))

  r_hlidar_twi_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlidar_twi.tif"))
  p_hlidar_twi_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlidar_twi.tif"))
  r_hlang_twi_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlang_twi.tif"))
  p_hlang_twi_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlang_twi.tif"))

  r_hlidar_tri_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlidar_tri.tif"))
  p_hlidar_tri_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlidar_tri.tif"))
  r_hlang_tri_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlang_tri.tif"))
  p_hlang_tri_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlang_tri.tif"))

  # ============================================================================
  # 2. DATA PREPARATION - HEIGHT AND TOPOGRAPHY
  # ============================================================================
  message("📊 Data preparation (H, Topo)...")

  tryCatch({
    # --- LiDAR height (reference) ---
    lidar <- rast(lidar_path) |> clamp(lower = 0, upper = 80, values = FALSE) # Clean
    vegh_lidar <- project(lidar, grid_30m, method = "average") |>             # Project and save
      save_raster(hlidar_file)

    # --- TWI (Topographic Wetness Index) ---
    twi_raw <- terra::crop(terra::rast(twi_30m_path), reg_extent) / 100
    twi <- mask(twi_raw, vegh_lidar) |>
      save_raster(twi_file)

    # --- 2.3 Lang et al. height ---
    vegh_lang_raw <- extent_to_tile_ids(reg_extent, tile_size = 3, return_raster = TRUE,
      source = "lang_vegh_10m", tiles_dir = vegh_10m_tiles_dir)

    vegh_lang <- raster_preprocess_save(
      input = vegh_lang_raw,
      output = hlang_file,
      target = vegh_lidar,
      na_value = 0,
      fun = mean,
      varname = "vegh_lang",
      if_aggregate = TRUE,
      if_resample = TRUE,
      if_mask = TRUE,
      if_return_raster = TRUE
    )

    # --- DEM, slope, aspect ---
    dem_raw <- extent_to_tile_ids( reg_extent, tile_size = 1, return_raster = TRUE,
      source = "COP30", tiles_dir = COP30_dir) |>
      resample(vegh_lidar) |>
      mask(vegh_lidar) |>
      save_raster(dem_file)

    topo_vars <- calc_tri(dem_raw, grid_30m, chunk_size, workers)

    # Save DEM, Slope, Aspect
    save_raster(topo_vars$dem, dem_file)
    save_raster(topo_vars$slope, slope_file)
    save_raster(topo_vars$aspect, aspect_file)
    save_raster(topo_vars$tri, tri_file)

    # ============================================================================
    # 4. CORRELATION ANALYSIS (500m windows)
    # ============================================================================

    message("\n📈 Calculating correlations...")
    stacked <- c(twi, tri, vegh_lang, vegh_lidar) # Stack all rasters
    df_win <- create_spatial_windows(
      stacked,
      value_vars = c("twi", "tri", "vegh_lang", "vegh_lidar"),
      dwin = 0.005
    ) # Create 500-m windows

    # ----------------------------------------------------------------------------
    # 4.1 Correlation with TWI
    # ----------------------------------------------------------------------------
    message("\n🔗 Correlation H-TWI...")
    df_cor_lang_twi <- calculate_correlation_bywin(df_win, x = "twi", y = "vegh_lang")
    df_cor_lidar_twi <- calculate_correlation_bywin(df_win, x = "twi", y = "vegh_lidar")

    df_combined_twi <- merge(
      df_cor_lang_twi, df_cor_lidar_twi,
      by = c("lon_mid", "lat_mid"),
      suffixes = c("_lang", "_lidar")
    )

    grid_500m <- create_aligned_template(vegh_lidar, res_out = 0.005, trim_input = TRUE)

    params_list <- list(
      var = c("correlation_lidar", "cor_pval_lidar", "correlation_lang", "cor_pval_lang"),
      file = c(r_hlidar_twi_file, p_hlidar_twi_file, r_hlang_twi_file, p_hlang_twi_file),
      name = c("r_hlidar_twi", "p_hlidar_twi", "r_hlang_twi", "p_hlang_twi")
    )

    mapply(function(var, file, name) {
      df_to_raster(
        df_combined_twi, "lon_mid", "lat_mid", var,
        template_raster = grid_500m,
        output_file = file,
        varname = name,
        overwrite = TRUE,
        return_raster = TRUE
      )
    }, var = params_list$var, file = params_list$file, name = params_list$name)

    # --- Calculate H-TRI correlation  ---
    message("\n🔗 Correlation H-TRI...")
    df_cor_lang_tri <- calculate_correlation_bywin(df_win, x = "tri", y = "vegh_lang")
    df_cor_lidar_tri <- calculate_correlation_bywin(df_win, x = "tri", y = "vegh_lidar")

    df_combined_tri <- merge(
      df_cor_lang_tri, df_cor_lidar_tri,
      by = c("lon_mid", "lat_mid"),
      suffixes = c("_lang", "_lidar")
    )

    # --- Save H-TRI correlation---
    params_list <- list(
      var = c("correlation_lidar", "cor_pval_lidar", "correlation_lang", "cor_pval_lang"),
      file = c(r_hlidar_tri_file, p_hlidar_tri_file, r_hlang_tri_file, p_hlang_tri_file),
      name = c("r_hlidar_tri", "p_hlidar_tri", "r_hlang_tri", "p_hlang_tri")
    )

    mapply(function(var, file, name) {
      df_to_raster(
        df_combined_tri, "lon_mid", "lat_mid", var,
        template_raster = r_temp_5km,
        output_file = file,
        varname = name,
        overwrite = TRUE,
        return_raster = TRUE
      )
    }, var = params_list$var, file = params_list$file, name = params_list$name)

    # ============================================================================
    # 5. SUMMARY STATISTICS
    # ============================================================================

    # --- 5.1 TWI correlation summary ---
    cat("\n=== H-TWI Summary ===\n")
    summary_twi <- data.frame(
      Variable = c("correlation_lang", "correlation_lidar"),
      Mean = c(mean(df_combined_twi$correlation_lang, na.rm = TRUE),
               mean(df_combined_twi$correlation_lidar, na.rm = TRUE)),
      SD = c(sd(df_combined_twi$correlation_lang, na.rm = TRUE),
             sd(df_combined_twi$correlation_lidar, na.rm = TRUE)),
      Min = c(min(df_combined_twi$correlation_lang, na.rm = TRUE),
              min(df_combined_twi$correlation_lidar, na.rm = TRUE)),
      Max = c(max(df_combined_twi$correlation_lang, na.rm = TRUE),
              max(df_combined_twi$correlation_lidar, na.rm = TRUE))
    )
    print(summary_twi)

    df_combined_twi$cor_diff <- df_combined_twi$correlation_lang - df_combined_twi$correlation_lidar
    cor_between_twi <- cor(
      df_combined_twi$correlation_lang,
      df_combined_twi$correlation_lidar,
      use = "complete.obs"
    )
    cat("Correlation between r(H_lang, TWI) and r(H_lidar, TWI):", cor_between_twi, "\n")

    # --- 5.2 TRI correlation summary ---
    cat("\n=== H-TRI Summary ===\n")
    summary_tri <- data.frame(
      Variable = c("correlation_lang", "correlation_lidar"),
      Mean = c(mean(df_combined_tri$correlation_lang, na.rm = TRUE),
               mean(df_combined_tri$correlation_lidar, na.rm = TRUE)),
      SD = c(sd(df_combined_tri$correlation_lang, na.rm = TRUE),
             sd(df_combined_tri$correlation_lidar, na.rm = TRUE)),
      Min = c(min(df_combined_tri$correlation_lang, na.rm = TRUE),
              min(df_combined_tri$correlation_lidar, na.rm = TRUE)),
      Max = c(max(df_combined_tri$correlation_lang, na.rm = TRUE),
              max(df_combined_tri$correlation_lidar, na.rm = TRUE))
    )
    print(summary_tri)

    df_combined_tri$cor_diff <- df_combined_tri$correlation_lang - df_combined_tri$correlation_lidar
    cor_between_tri <- cor(
      df_combined_tri$correlation_lang,
      df_combined_tri$correlation_lidar,
      use = "complete.obs"
    )
    cat("Correlation between r(H_lang, tri) and r(H_lidar, tri):", cor_between_tri, "\n")

    message("\n🎉 Validation completed successfully for: ", file_base)
    return(TRUE)

  }, error = function(e) {
    message("\n❌ Error in validation for ", file_base, ": ", e$message)
    return(FALSE)
  })
}

# ----------------- Parallel execution for all tiles-----------------
# Set up cluster plan
plan(cluster, workers = workers)

# Run in parallel
tictoc::tic("🚀 Parallel processing of tiles")

results <- future_map(
  file_bases,
  validate_vegetation_heights,
  .progress = FALSE,
  .options = furrr_options(seed=TRUE)
)

plan(sequential)
tictoc::toc()

# Summarize results
success_count <- sum(unlist(results))
fail_count <- length(results) - success_count
message(sprintf("✅ Completed: %d succeeded, ❌ %d failed.", success_count, fail_count))

# ------ Single tile check (optional) ---------------------------------------
# validate_vegetation_heights(file_bases[1])
# validate_vegetation_heights(file_bases[11])
