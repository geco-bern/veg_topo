
# ==============================================================================
# Moisture Index Aggregation
# Purpose:
#   Aggregate global 0.05° (5km) moisture index data to 0.5° (55km) resolution.
#
# Run time:
# ~ 1 min on UBELIX
# ==============================================================================

# -------------------- 1. Setup Environment ------------------------------------
library(terra)

source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/create_aligned_template.R"))

# Create output directory
if (!dir.exists(dirname(mi_55km_file))) dir.create(dirname(mi_55km_file), recursive = TRUE)

# ---------------- 2. Zonal aggregation (0.05° → 0.5°) -----------------
message("Aggregation MI (0.05° → 0.5°)")

# Create template raster aligned to 55km grid and zonal aggregate
align_template_55km <- create_aligned_template(twi_450m_mosaic_clean_path, res_out = 0.5)
raster_preprocess_save(
  input   = mi_5km_file,
  output  = mi_55km_file,
  target  = align_template_55km,
  varname = "mi",
  if_zonal = TRUE,
  fun = "mean",
  if_aggregate = FALSE,
  if_resample    = FALSE,
  if_return_raster = FALSE
)

# # -------------------- 4. Quick Check & Visualization --------------------------

# mi_55km_r <- terra::rast(mi_55km_file)
# print(mi_55km_r)
# summary(mi_55km_r)
# plot(mi_55km_r, main = "Moisture Index (0.5°)")

# mi_5km_r  <- terra::rast(mi_5km_file)
# print(mi_5km_r)
# summary(mi_5km_r)
# plot(mi_5km_r, main = "Moisture Index (0.05°)")





