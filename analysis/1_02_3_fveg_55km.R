# ==============================================================================
# Vegetated Area Fraction Aggregation (450m → 0.5°)
#
# Purpose:
#   Aggregate high-resolution fraction of vegetated area (450m) data to 0.5° (~55km)
#   resolution by calculating the fraction of vegetated area within each
#   coarse grid cell (non-zero vegetation height / total pixels).
# ==============================================================================

# -------------------- 1. Setup Environment ------------------------------------
library(terra)
library(exactextractr) # Zonal aggregation

# other custom functions
source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/create_aligned_template.R"))

# Create output directory
if (!dir.exists(dirname(fveg_55km_path))) dir.create(dirname(fveg_55km_path), recursive = TRUE)

# -------------------- 2. Aggregation & Resampling -----------------------------
message("Aggregating vegetated area fraction...")

# Create template raster aligned to 5km grid
align_template_55km <- create_aligned_template(twi_450m_mosaic_clean_path, res_out=0.5)

# real vegetated fraction based on 10m height data
raster_preprocess_save(
  input        = fveg_real_450m_mosaic_path, # calculate from 10-m height data
  output       = fveg_real_55km_path,
  target       = align_template_55km,
  varname      = "fveg",
  if_zonal = TRUE,
  fun = function(values, coverage_fractions) {
    if (all(is.na(values))) return(NA_real_)
    # Area-weighted mean (more accurate for partial pixels)
    sum(values * coverage_fractions, na.rm = TRUE) / sum(coverage_fractions, na.rm = TRUE)
  },
  if_aggregate = FALSE,
  if_resample  = FALSE,
  if_return_raster = FALSE
)

# vegetated area fraction based on valid height data on vegh_450m_mosaic_path
# The weighted proportion of a target pixel's area that is covered by non-NA (valid) source pixels.
raster_preprocess_save(
  input        = vegh_450m_mosaic_path,
  output       = fveg_55km_path,
  target       = align_template_55km,
  varname      = "fveg",
  if_zonal = TRUE,
  fun = function(values, coverage_fractions) {
    if (all(is.na(values))) return(NA_real_)
    # Count of non-NA values weighted by coverage fractions
    sum((!is.na(values)) * coverage_fractions, na.rm = TRUE) / sum(coverage_fractions, na.rm = TRUE)
  },
  if_aggregate = FALSE,
  if_resample  = FALSE,
  if_return_raster = FALSE
)

# # -------------------- 3. Optional Check ---------------------------------------
# # Uncomment this section to verify results visually
# message("Checking input and output rasters...")
# r_in  <- rast(fveg_real_450m_mosaic_path)
# r_out <- rast(fveg_55km_path)
# r_out2 <- rast(fveg_real_55km_path)
#
# print(r_in)
# print(r_out)
#
# par(mfrow = c(1, 2))
# plot(r_in,  main = "Vegetation Height (450m)")
# plot(r_out, main = "Vegetated Area Fraction (0.5°)")
# plot(r_out2, main = "Vegetation Cover Fraction (0.5°)")
#
# par(mfrow = c(1, 1))
#
# summary(r_out)

