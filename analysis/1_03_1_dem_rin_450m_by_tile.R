# ~ 15.5 h: 50 cores, 300 G memory on UBELIX
# =============== Setup and configuration ======================================
library(terra)
library(dplyr)
library(tidyr)
library(fs)
library(stringr)
library(parallel)
library(future)
library(furrr)
library(meteoland)

source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/df_to_raster.R"))
source(here::here("R/create_aligned_template.R"))
source(here::here("R/cacl_meteoland_sw_in.R"))
source(here::here("R/convert_lat.R"))
source(here::here("R/convert_lon.R"))
source(here::here("R/merge_dem_neighbors.R"))
source(here::here("R/calc_tri.R"))

# --- PARALLEL CONFIGURATION (by hostname) ---
# Set worker configuration based on system
if (hostname == "dash") {
  chunk_size <- 1000
  workers <- 16
} else {
  chunk_size <- 10000
  workers <- 100
}
message("→ Using ", workers, " workers and chunk_size = ", chunk_size, " for parallel processing.\n")

# --- Total tiles: 19429 ---
# Test mode (recommended for first run)
start_idx <- 1
end_idx   <- 19429

# Full run (uncomment when ready)
# start_idx <- 1
# end_idx   <- 19429

# --- Input and output directory ---
if (!dir.exists(COP30_dir)) stop(paste("Input directory does not exist:", COP30_dir))

if (!dir.exists(rin_450m_tiles_dir)) dir.create(rin_450m_tiles_dir, recursive = TRUE)
if (!dir.exists(dem_450m_tiles_dir)) dir.create(dem_450m_tiles_dir, recursive = TRUE)
if (!dir.exists(slope_450m_tiles_dir)) dir.create(slope_450m_tiles_dir, recursive = TRUE)
if (!dir.exists(aspect_450m_tiles_dir)) dir.create(aspect_450m_tiles_dir, recursive = TRUE)

# =============== Single tile processing function ==============================

cal_rin_by_tile <- function(lat, lon) {

  # Output file paths per tile
  reg_id <- paste0(convert_lat(lat), "_", convert_lon(lon))

  rin_file <- file.path(rin_450m_tiles_dir,
                            paste0("radiation_index_", reg_id, "_15_arcscd.tif"))
  dem_file <- file.path(dem_450m_tiles_dir,
                            paste0("dem_", reg_id, "_15_arcscd.tif"))

  slope_file <- file.path(slope_450m_tiles_dir,
                            paste0("slope_", reg_id, "_15_arcscd.tif"))

  aspect_file <- file.path(aspect_450m_tiles_dir,
                            paste0("aspect_", reg_id, "_15_arcscd.tif"))

  # Skip if already processed
  if (file.exists(rin_file) && file.exists(dem_file)) return(TRUE)

  tryCatch({

    # Create alignment grid (30-m /1-arcsec resolution, )
    grid_30m <- create_aligned_template(ext(lon, lon + 1, lat, lat + 1), res_out = 1/3600)

    # --- Data pre: dem, slope and aspect ---

    # Load DEM with neighbors
    dem_ex <- merge_dem_neighbors(lat, lon, file_dir = COP30_dir)
    if (is.null(dem_ex)) stop("Failed to load DEM data for this tile")

    # Compute slope (450m resolution)
    slope_450m <- terrain(dem_ex, v = "slope", unit = "degrees") |>
      resample(grid_30m, method = "bilinear") |>
      aggregate(c(15, 15)) |>
      save_raster(slope_file)

    # Compute aspect (450m resolution, circular mean)
    aspect <- terrain(dem_ex, v = "aspect", unit = "radians")
    aspect_cos <- cos(aspect) |>
      resample(grid_30m, method = "bilinear") |>
      aggregate(c(15, 15))
    aspect_sin <- sin(aspect) |>
      resample(grid_30m, method = "bilinear") |>
      aggregate(c(15, 15))
    aspect_450m <- (atan2(aspect_sin, aspect_cos) * 180 / pi) %% 360 |>
      save_raster(aspect_file)

    # Resample and save DEM
    dem_450m <- resample(dem_ex, slope_450m) |>
      save_raster(dem_file)

    # Extract topography to dataframe
    df_topo <- as.data.frame(slope_450m, xy = TRUE) |>
      left_join(as.data.frame(aspect_450m, xy = TRUE), by = c("x", "y")) |>
      tidyr::drop_na()
    colnames(df_topo) <- c("lon", "lat", "slope", "aspect")

    # --- Radiation calculation ---
    # Slope surface radiation
    sw_meteoland_surf_vec <- cacl_meteoland_sw_in(df_topo$lat, df_topo$slope, df_topo$aspect, 2020)

    # Flat surface radiation (cached by latitude)
    sw_meteoland_flat_vec <- unlist(ave(df_topo$lat, df_topo$lat,
                                        FUN = function(x) cacl_meteoland_sw_in(x[1], 0, 0, 2020)))

    # Calculate and save Topographic Radiation Index
    rin_meteoland_vec <- sw_meteoland_surf_vec / sw_meteoland_flat_vec
    df_calc <- df_topo[, c("lon", "lat")] |>
      mutate(rin = rin_meteoland_vec)

    rin_450m <- suppressMessages(
      df_to_raster(df_calc, "lon", "lat", "rin", slope_450m, output_file = rin_file)
    )

    return(TRUE)

  }, error = function(e) {
    message(sprintf("❌ Tile (%s, %s) failed: %s", convert_lat(lat), convert_lon(lon), e$message))
    return(FALSE)
  })
}

# =============== Load raw data (dem) ==========================================
dem_files_all <- fs::dir_ls(
  path = COP30_dir,
  glob = "*_DEM.tif",
  recurse = TRUE
)
message(sprintf("Found %d DEM tiles", length(dem_files_all)))

# Extract lat/lon from filenames
lat_str <- substring(str_extract(dem_files_all, "_[NS]\\d{2}"), 2)
lon_str <- substring(str_extract(dem_files_all, "_[EW]\\d{3}"), 2)
lats <- unname(sapply(lat_str, convert_lat))
lons <- unname(sapply(lon_str, convert_lon))

# Filter high-latitude tiles
keep_idx <- lats >= -57 & lats <= 87
dem_files_all <- dem_files_all[keep_idx]
lats <- lats[keep_idx]
lons <- lons[keep_idx]

message(sprintf("After filtering high-Lat tiles,\n  Found %d DEM tiles", length(dem_files_all)))

# --- Start process ---
message(sprintf("⭐️⭐️⭐️ Processing tiles %d to %d out of %d total ⭐️⭐️⭐️", start_idx, end_idx, length(dem_files_all)))

# =============== Parallel execution ===========================================

plan(cluster, workers = workers)

# Get tile indices
tile_indices <- start_idx:end_idx
n_tiles <- length(tile_indices)

message("========================================")
message(sprintf("HPC Mode (Silent):"))
message(sprintf("  Total tiles: %d", n_tiles))
message(sprintf("  Workers: %d", workers))
message("========================================")

tictoc::tic("🚀 HPC Parallel processing")

results <- future_map(tile_indices, function(idx) {
  cal_rin_by_tile(lats[idx], lons[idx])
}, .progress = FALSE, .options = furrr_options(seed = TRUE))

plan(sequential)
tictoc::toc()

# --- summary message ---
all_results <- unlist(results)
success_count <- sum(all_results, na.rm = TRUE)
fail_count <- length(all_results) - success_count

message("========================================")
message(sprintf("✅ %d succeeded, ❌ %d failed (%.1f%%)",
                success_count, fail_count, success_count / n_tiles * 100))
message("========================================")

# # =============== One tile for step by step test =============================
# id_test <- 22
# lat <- lats[[id_test]]
# lon <- lons[[id_test]]
#
# cal_rin_by_tile(lat, lon)
#
# plot(rast("/data_2/scratch/ting/veg_topo_data/data/global_dem_slope_aspect_450m/dem_1_1_deg/dem_N00_E029_15_arcscd.tif"))
# plot(rast("/data_2/scratch/ting/veg_topo_data/data/global_rin_450m/1_1_deg/radiation_index_N00_E029_15_arcscd.tif"))
#
# # ---- 单片测试参数设置 ----
# # 选择一个测试瓦片
# id_test <- 11111
# lat <- lats[[id_test]]
# lon <- lons[[id_test]]
#
# # 输出文件路径
# rin_file <- file.path(rin_450m_tiles_dir,
#                           paste0("radiation_index_", convert_lat(lat), "_", convert_lon(lon), "_15_arcscd.tif"))
# dem_file <- file.path(dem_450m_tiles_dir,
#                           paste0("dem_", convert_lat(lat), "_", convert_lon(lon), "_15_arcscd.tif"))
#
# # 创建对齐网格
# grid_30m <- create_aligned_template(ext(lon, lon + 1, lat, lat + 1), res_out = 1/3600)
#
# # ---- 处理 DEM ----
# # 合并相邻 DEM 瓦片
# dem_ex <- merge_dem_neighbors(lat, lon, file_dir = COP30_dir)
#
# # 检查 DEM 是否有效
# if (is.null(dem_ex)) {
#   stop("Failed to load DEM data for this tile")
# }


