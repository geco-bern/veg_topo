# ==============================================================================
# Ecoregion Vector to Raster Pipeline with Area-Weighted Statistics
#
# Workflow:
#   1. Create high-resolution template (0.01° for accurate area calc)
#   2. Rasterize ecoregions at high resolution
#   3. Aggregate to target resolution (0.05°) with area-weighted mean
#   4. Output: proportion of each biome per pixel (continuous values)
#
# Runtime:
# ~43 mins on UBELIX
# ==============================================================================

# -------------------- 1. Setup ------------------------------------------------
library(terra)
library(exactextractr)
library(sf)

source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/create_aligned_template.R"))

# Create output directory
if (!dir.exists(dirname(ecoregion_5km_path))) dir.create(dirname(ecoregion_5km_path), recursive = TRUE)

# -------------------- 2. High-resolution rasterization ------------------------
message("Creating high-resolution ecoregion raster (0.01° for accuracy)...")

# Create 0.05° template aligned with existing raster
r_tar <- create_aligned_template(twi_450m_mosaic_clean_path)

# Create high-res template (5x finer than target)
r_highres_template <- rast(ext(r_tar), resolution = 0.01, crs = crs(r_tar))

# Load ecoregion vector
v <- vect(ecoregion_path)

# Rasterize at high resolution (preserve categorical codes)
r_highres <- rasterize(v, r_highres_template, field = "BIOME_NUM")

# -------------------- 3. Aggregate with area weighting ------------------------
message("Aggregating to target resolution (area-weighted proportion)...")
r_out <- raster_preprocess_save(
  input        = r_highres,
  output       = ecoregion_5km_path,
  target       = r_tar,
  varname      = "fpa",
  if_zonal     = TRUE,
  fun          = "majority",
  if_aggregate = FALSE,
  if_resample  = FALSE,
  if_return_raster = FALSE
)
