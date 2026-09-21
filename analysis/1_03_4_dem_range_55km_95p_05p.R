# ============================================================
# Elevation Range Calculation (0.5° grid from 450m DEM)
#
# Steps:
#   1. Load DEM and target grid
#   2. Compute 95th and 5th percentile elevation per grid cell
#   3. Compute elevation range (95p - 05p)

# about ~ 11 min on UBELIX
# ============================================================

library(terra)
library(sf)
source(here::here("R/config.R"))
source(here::here("R/create_aligned_template.R"))
source(here::here("R/raster_preprocess_save.R"))

# Create output directory
if (!dir.exists(dirname(dem_rg_95p_05p_55km_path))) dir.create(dirname(dem_rg_95p_05p_55km_path), recursive = TRUE)

# ---------------- 1. Create template grid ------------------------
# Create template raster aligned to 55km grid
align_template_55km <- create_aligned_template(twi_450m_mosaic_clean_path, res_out = 0.5)

# ---------------- 2. Compute 95th percentile -----------------
message("Calculating 95th percentile elevation...")
dem_95p <- raster_preprocess_save(
  input   = dem_450m_mosaic_path,
  output  = NULL,
  target  = align_template_55km,
  varname = "dem_95p",
  if_zonal = TRUE,
  if_aggregate = FALSE,
  fun = function(values, coverage_fractions) {
    # Keep only pixels with at least 10% coverage to reduce edge noise
    keep <- coverage_fractions >= 0.1
    values_filtered <- values[keep]
    # Return NA if insufficient pixels for reliable quantile
    if (length(values_filtered) < 2) {
      return(NA_real_)
    }
    as.numeric(quantile(values_filtered, 0.95, na.rm = TRUE))
  },
  if_resample    = FALSE,
  if_return_raster = TRUE
)

# ---------------- 3. Compute 5th percentile ------------------
message("Calculating 05th percentile elevation...")
dem_05p <- raster_preprocess_save(
  input   = dem_450m_mosaic_path,
  output  = NULL,
  target  = align_template_55km,
  varname = "dem_05p",
  if_zonal = TRUE,
  fun = function(values, coverage_fractions) {
    # Keep only pixels with at least 10% coverage to reduce edge noise
    keep <- coverage_fractions >= 0.1
    values_filtered <- values[keep]
    # Return NA if insufficient pixels for reliable quantile
    if (length(values_filtered) < 2) {
      return(NA_real_)
    }
    as.numeric(quantile(values_filtered, 0.05, na.rm = TRUE))
  },
  if_aggregate = FALSE,
  if_round_fact  = TRUE,
  if_resample    = FALSE,
  if_return_raster = TRUE
)

# ---------------- 4. Compute elevation range -----------------
message("Computing elevation range (95p - 05p)...")
dem_rg <- dem_95p - dem_05p

terra::writeCDF(dem_rg, dem_rg_95p_05p_55km_path, overwrite = TRUE, varname = "dem_rg_95p_05p")
if (file.exists(dem_rg_95p_05p_55km_path)) message("✅ Saved : ", dem_rg_95p_05p_55km_path)

# # ---------------- 5 (optional). Check the output ----------
# r <- terra::rast(dem_rg_95p_05p_55km_path)
# r
# summary(r)
# plot(r, main = "Elevation range (0.5°)")
