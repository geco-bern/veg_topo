library(terra)
library(dplyr)
library(sf)
library(tidyr)
library(furrr)
library(future)
library(ggplot2)
library(patchwork)

source(here::here("R/config.R"))
source(here::here("R/get_lonlat_extent.R")) # in create_aligned_template() to process lidar raster
source(here::here("R/create_aligned_template.R"))
source(here::here("R/create_spatial_windows.R"))
source(here::here("R/calculate_correlation_bywin.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/df_to_raster.R"))
source(here::here("R/aggregate_topography.R"))
source(here::here("R/extent_to_tile_ids.R"))
source(here::here("R/helpers.R"))
source(here::here("R/calc_meteoland_sw_in.R"))

# Set worker numbers for different system
hostname <- trimws(tolower(system("hostname", intern = TRUE)))
if (hostname == "dash") workers <- 4 else workers <- 12
message("→ using ", workers, " workers")

# All ALS_MAX data
file_bases <- c("30SWH","30TWN","30UVD","31UFT","32TMS",
                "32TMT", "34WFS", "35VMF","35VNL","35WMR","35WNT",
                "06VXR", "08WNA", "12VUN", "16PHS", "32MPC", "32MQE",
                "06WVT", "10TER", "13UEA", "16SGE", "32MPE", "32NNF")

# get target resolution
twi <- rast(twi_450m_mosaic_clean_path)
res_tar <- res(twi)

# Function to process for each Lidar tile
validate_vegetation_heights <- function(file_base) {

  # ============================================================================
  # Initialization
  # ============================================================================
  message("⭐️⭐️⭐️ Processing: ", file_base, " ⭐️⭐️⭐️")

  # Create output directory
  if (!dir.exists(h_validation_dir)) dir.create(h_validation_dir, recursive = TRUE)

  # ----- Check if input file exists -----
  # LiDAR input path and alignment grid
  lidar_path <- file.path(lidar_asl_dir, paste0(file_base, ".tif"))
  lidar_path <- if (file.exists(lidar_path)) lidar_path else file.path(lidar_lvis_dir, paste0(file_base, ".tif"))
  if (!file.exists(lidar_path)) return(message("❌ LiDAR file not found"))

  # ----- Region info -----
  # Get the extent using alignment template
  align_r <- create_aligned_template(lidar_path, res_out = res_tar[1], snap_step = 0.025 )
  reg_extent <- ext(align_r)
  message("Extent:",reg_extent[1], ", ", reg_extent[2],", ",reg_extent[3],", ",reg_extent[4], " (xmin,xmax,ymin,ymax)")

  # ----- Output file paths -----
  hlidar_file <- file.path(h_validation_dir, paste0(file_base, "_hlidar_450m.tif"))
  hlang_file <- file.path(h_validation_dir, paste0(file_base, "_hlang_450m.tif"))
  twi_file <- file.path(h_validation_dir, paste0(file_base, "_twi_450m.tif"))
  dem_file <- file.path(h_validation_dir, paste0(file_base, "_dem_450m.tif"))
  slope_file <- file.path(h_validation_dir, paste0(file_base, "_slope_450m.tif"))
  aspect_file <- file.path(h_validation_dir, paste0(file_base, "_aspect_450m.tif"))
  tri_file <- file.path(h_validation_dir, paste0(file_base, "_tri_450m.tif"))

  r_hlidar_twi_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlidar_twi_5km.tif"))
  p_hlidar_twi_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlidar_twi_5km.tif"))
  r_hlang_twi_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlang_twi_5km.tif"))
  p_hlang_twi_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlang_twi_5km.tif"))

  r_hlidar_tri_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlidar_tri_5km.tif"))
  p_hlidar_tri_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlidar_tri_5km.tif"))
  r_hlang_tri_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlang_tri_5km.tif"))
  p_hlang_tri_file <- file.path(h_validation_dir, paste0(file_base, "_p_hlang_tri_5km.tif"))

  # ============================================================================
  # Data pre: height and topo
  # ============================================================================
  message("📊 Data preparation (H, Topo)...")

  tryCatch({
    # --- LiDAR height (reference) ---
    # load, clean, project and save
    lidar <- rast(lidar_path) |> clamp(lower = 0, upper = 70, values = FALSE) # Clean
    vegh_lidar <- project(lidar, align_r, method = "average") |>
      save_raster(hlidar_file)

    # --- TWI (Topographic Wetness Index) ---
    twi <- terra::crop(terra::rast(twi_450m_mosaic_clean_path), reg_extent) |>
      mask(vegh_lidar) |>
      save_raster(twi_file)

    # --- height from Lang et al. 2023 ---
    vegh_lang <- terra::crop(terra::rast(vegh_450m_mosaic_path), reg_extent) |>
      mask(vegh_lidar) |>
      save_raster(hlang_file)

    # --- DEM, slope and aspect ---
    # crop and save
    dem <- terra::crop(terra::rast(dem_450m_mosaic_path), reg_extent) |>
      save_raster(dem_file)
    slope <- terra::crop(terra::rast(slope_450m_mosaic_path), reg_extent) |>
      save_raster(slope_file)
    aspect <- terra::crop(terra::rast(aspect_450m_mosaic_path), reg_extent) |>
      save_raster(aspect_file)
    rin <- terra::crop(terra::rast(tri_450m_path), reg_extent) |>
      save_raster(tri_file)

    # ============================================================================
    # Pearson's correlation ( 5-km fixed window)
    # ============================================================================
    # --- Date pre and fixed windows ---
    message("\n📈 Calculating correlations...")
    stacked <- c(twi, rin, vegh_lang, vegh_lidar)

    message("\n🪟 Create spatial windows...")
    # Create spatial windows
    df_win <- create_spatial_windows(
      stacked,
      value_vars = c("twi", "tri", "vegh_lang", "vegh_lidar"),
      dwin = 0.05
    )

    # --- Calculate H-TWI correlation  ---
    message("\n🔗 Correlation H-TWI...")
    df_cor_lang_twi <- calculate_correlation_bywin(df_win, x = "twi", y = "vegh_lang")
    df_cor_lidar_twi <- calculate_correlation_bywin(df_win, x = "twi", y = "vegh_lidar")

    df_combined_twi <- merge(
      df_cor_lang_twi, df_cor_lidar_twi,
      by = c("lon_mid", "lat_mid"),
      suffixes = c("_lang", "_lidar")
    )

    # --- Save H-TWI correlation---
    # aligned grid
    r_temp_5km <- create_aligned_template(vegh_lidar, res_out = 0.05)

    params_list <- list(
      var = c("correlation_lidar", "cor_pval_lidar", "correlation_lang", "cor_pval_lang"),
      file = c(r_hlidar_twi_file, p_hlidar_twi_file, r_hlang_twi_file, p_hlang_twi_file),
      name = c("r_hlidar_twi", "p_hlidar_twi", "r_hlang_twi", "p_hlang_twi")
    )

    mapply(function(var, file, name) {
      df_to_raster(
        df_combined_twi, "lon_mid", "lat_mid", var,
        template_raster = r_temp_5km,
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

    # --- 5.2 Rin correlation summary ---
    cat("\n=== H-Rin Summary ===\n")
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
    cat("Correlation between r(H_lang, Rin) and r(H_lidar, Rin):", cor_between_tri, "\n")

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
