# ==============================================================================
# Protected Areas (PA) Rasterization and Aggregation (0.02° → 0.5°)
#
# Workflow:
#   1. Rasterize multiple PA shapefiles to 0.02° grid (binary: 1=protected)
#   2. Merge rasters using maximum operator (union of all PAs)
#   3. Aggregate to 0.5° using area-weighted mean (fractional coverage)
#
# Runtime:
# ~5 mins on UBELIX
# ==============================================================================

# -------------------- 1. Setup ------------------------------------------------
library(terra)
library(sf)

source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/create_aligned_template.R"))

# Create output directory
if (!dir.exists(dirname(fpa_55km_path))) dir.create(dirname(fpa_55km_path), recursive = TRUE)

message("Starting PA rasterization and aggregation...")

# -------------------- 2. File paths -------------------------------------------
message("Configuring file paths...")

# Input shapefiles
shp_files <- c(pa_shp0, pa_shp1, pa_shp2)

# Intermediate files
temp_dir <- tempdir()
rasterized_files <- file.path(temp_dir, paste0("pa_", 0:2, ".tif"))
merged_22km_path <- file.path(temp_dir, "pa_merged_22km.tif")

# -------------------- 3. Rasterize shapefiles to 0.02° ------------------------
message("Rasterizing shapefiles to 0.02° grid...")

# -burn 1 -> 1 as protected
rasterize_shp <- function(shp_path, out_tif) {
  cmd <- sprintf(
    "gdal_rasterize -burn 1 -tr 0.02 0.02 -te -180 -60 180 90 -ot Byte -of GTiff '%s' '%s'",
    shp_path, out_tif
  )
  system(cmd)
}

mapply(rasterize_shp, shp_files, rasterized_files)

# -------------------- 4. Merge rasters (union) --------------------------------
message("Merging rasters with maximum operator (union of all PAs)...")

r1 <- rast(rasterized_files[1])
r2 <- rast(rasterized_files[2])
r3 <- rast(rasterized_files[3])

# merge by terra::max
r_merged_terra <- max(r1, r2, r3, na.rm = TRUE)

# -------------------- 5. Aggregate to 0.5° (area-weighted fraction) -----------
message("Aggregating to 0.5° grid (fractional PA coverage)...")

# Create 0.5° template aligned with existing raster
align_template_55km <- create_aligned_template(twi_450m_mosaic_clean_path, res_out = 0.5)

# Area-weighted aggregation: accounts for partial pixel coverage at boundaries
# Why 0 treat as NA here
r_out <- raster_preprocess_save(
  input        = r_merged_terra,
  output       = fpa_55km_path,
  target       = align_template_55km,
  varname      = "fpa",
  if_zonal     = TRUE,
  fun = function(values, coverage_fractions) {
    if (all(is.na(values))) return(0)
    # Area-weighted mean (more accurate for partial pixels)
    sum(values * coverage_fractions, na.rm = TRUE) / sum(coverage_fractions, na.rm = TRUE)
  },
  if_aggregate = FALSE,
  if_resample  = FALSE,
  if_return_raster = FALSE
)

# # -------------------- 6. Optional inspection ---------------------------------
# r_check <- rast(fpa_55km_path)
# print(r_check)
# summary(r_check)
# plot(r_check, main = "Protected Area Fraction (0.5°)")
