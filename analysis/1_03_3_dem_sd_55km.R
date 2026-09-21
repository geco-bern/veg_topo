# ==============================================================================
# DEM Standard Deviation Aggregation (450m → 0.5°)
#
# Purpose:
#   Aggregate high-resolution DEM (450m) data to 0.5° (55km) resolution by
#   calculating standard deviation (sd) within each coarse grid cell.
# ==============================================================================

# -------------------- 1. Setup Environment ------------------------------------
library(terra)

source(here::here("R/config.R"))
source(here::here("R/create_aligned_template.R"))
source(here::here("R/raster_preprocess_save.R"))

# Create output directory
if (!dir.exists(dirname(dem_sd_55km_path))) dir.create(dirname(dem_sd_55km_path), recursive = TRUE)

# -------------------- 2. Aggregation & Resampling -----------------------------
message("Aggregating DEM using standard deviation (sd)...")

# Create template raster aligned to 55km grid
align_template_55km <- create_aligned_template(twi_450m_mosaic_clean_path, res_out = 0.5)

raster_preprocess_save(
  input        = dem_450m_mosaic_path,
  output       = dem_sd_55km_path,
  target       = align_template_55km,
  varname      = "dem_sd",
  if_zonal = TRUE,
  fun          = "stdev",
  if_aggregate = FALSE,
  if_resample  = FALSE,
  if_return_raster = FALSE
)

message("✅ Completed.")

# # -------------------- 3. Optional Check ---------------------------------------
# # Uncomment this section to verify results visually
# message("Checking input and output rasters...")
# r_in  <- rast(dem_450m_mosaic_path)
# r_out <- rast(dem_sd_55km_path)
#
# print(r_in)
# print(r_out)
# summary(r_out)
#
# par(mfrow = c(1, 2))
# plot(r_in,  main = "DEM (450m)")
# plot(r_out, main = "DEM SD (0.5°)")
# par(mfrow = c(1, 1))


message("✅ Visualization complete.")
