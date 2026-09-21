# ==============================================================================
# Purpose:
# Compute mean annual climatic variables from monthly WorldClim .tif files (1km)
# and zonal aggregate to 5km (0.05°):
#   - Mean Annual Temperature (MAT)  - monthly mean
#   - Annual Total Precipitation (MAP) - monthly sum
#   - Annual Total Solar Radiation (SRAD) - monthly sum

# Run time:
# ~132 mins on UBELIX
# ==============================================================================

# ------------- Setup ----------------------------------------------------------
library(terra)
library(sf)
source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/create_aligned_template.R"))

# Create output directory
if (!dir.exists(dirname(mat_5km_file))) dir.create(dirname(mat_5km_file), recursive = TRUE)

# ------------- Process: MAT (mean), MAP & SRAD (sum) --------------------------
vars    <- c("tavg", "prec", "srad")
outputs_1km <- c(mat_1km_file, map_1km_file, srad_1km_file)
# outputs_5km <- c(mat_5km_file, map_5km_file, srad_5km_file)
# outputs_450m <- c(mat_450m_file, map_450m_file, srad_450m_file)

# Create 5km template once (reused for all variables)
align_template_5km <- create_aligned_template(twi_450m_mosaic_clean_path)
for (i in seq_along(vars)) {
  var    <- vars[i]
  message("⭐️⭐️⭐️ Processing: ", var, " ⭐️⭐️⭐️")

  output_1km <- outputs_1km[i]
  # output_5km <- outputs_5km[i]
  # output_450m <- outputs_450m[i]

  # more info here:
  # https://api.rdocumentation.org/packages/geodata/versions/0.5-8/topics/worldclim

  # Monthly file paths: wc2.1_30s_{var}_01.tif to ..._12.tif
  # Load 12 months and compute annual statistic
  monthly_rasters <- rast(paste0(worldclim_1km_dir, "wc2.1_30s_", var, "_", sprintf("%02d", 1:12), ".tif"))
  annual_result <- switch(var,
                          tavg = mean(monthly_rasters),  # °C, annual mean temperature
                          prec = sum(monthly_rasters),   # mm, annual total precipitation
                          srad = mean(monthly_rasters)   # kJ m⁻² day⁻¹, mean monthly daily radiation or annual mean daily solar radiation
  )

  # # Save to temporary file
  message("Saving 1-km...")
  raster_preprocess_save(
    input            = annual_result,
    output           = output_1km,
    target           = NULL,
    varname          = var,
    if_zonal         = FALSE,
    if_aggregate     = FALSE,
    if_resample      = FALSE,
    if_return_raster = FALSE
  )

  # # Align to 5km grid and aggregate zones
  # message("Resampling to 5-km...")
  # raster_preprocess_save(
  #   input            = annual_result,
  #   output           = output_5km,
  #   target           = align_template_5km,
  #   varname          = var,
  #   if_zonal         = TRUE,
  #   fun              = "mean", # for exact_extract()
  #   if_aggregate     = FALSE,
  #   if_resample      = FALSE,
  #   if_return_raster = FALSE
  # )
  #
  # # Resample to 450m grid (as Marthews et al.'s TWI data)
  # message("Resampling to 450-m...")
  # raster_preprocess_save(
  #   input            = annual_result,
  #   output           = output_450m,
  #   target           = twi_450m_mosaic_clean_path,
  #   varname          = var,
  #   if_zonal         = FALSE,
  #   if_aggregate     = FALSE,
  #   if_resample      = TRUE,
  #   if_return_raster = FALSE
  # )
}

# -------------------- Optional Check ---------------------------
# r <- rast(mat_5km_file)
# print(r)
# summary(r)
# plot(r, main = "Mean annual temperature (5km)")
#
# r <- rast(map_5km_file)
# print(r)
# summary(r)
# plot(r, main = "Mean annual precipitation (5km)")
#
# r <- rast(srad_5km_file)
# print(r)
# summary(r)
# plot(r, main = "Mean annual solar radiation (5km)")
