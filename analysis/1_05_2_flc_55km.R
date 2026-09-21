# ==============================================================================
# Fused Data Aggregation and Resampling (0.05° → 0.5°)
#
# Purpose:
#   Aggregate and resample fused 5km (0.05°) raster data to 0.5° (55km) resolution.
#
# Run time:
# ~ 2 mins on UBELIX
# ==============================================================================

# -------------------- 1. Setup Environment ------------------------------------
library(terra)

source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/create_aligned_template.R"))

# Create output directory
if (!dir.exists(dirname(fused_55km_file))) dir.create(dirname(fused_55km_file), recursive = TRUE)

# ---------------- 2. Zonal aggregation (0.05° → 0.5°) -----------------
message("Aggregation (0.05° → 0.5°)")

# Set all files
input_files <- c(fused_5km_file, fbare_5km_file, fwater_5km_file, fsnow_5km_file)
output_files <- c(fused_55km_file, fbare_55km_file, fwater_55km_file, fsnow_55km_file)
var_names <- c("fused", "fbare", "fwater", "fsnow")

# Create template raster aligned to 55km grid and zonal aggregate
align_template_55km <- create_aligned_template(twi_450m_mosaic_clean_path, res_out = 0.5)
mapply(function(input, output, varname) {
  raster_preprocess_save(
    input   = input,
    output  = output,
    target  = align_template_55km,
    varname = varname,
    if_zonal = TRUE,
    fun = "mean",
    if_aggregate = FALSE,
    if_resample = FALSE,
    if_return_raster = FALSE
  )
}, input_files, output_files, var_names)

# # -------------------- 4. Quick Check & Visualization --------------------------
# fused_55km_r <- terra::rast(fused_55km_file)
# print(fused_55km_r)
# summary(fused_55km_r)
# plot(fused_55km_r, main = "Used land fraction (0.5°)")

