# ------Set up------------------------------------------------------------------
library(terra)
source(here::here("R/config.R"))

# ------Data Clean: 10 min--------------------------------------------------------------

# Load the raster
twi_r <- terra::rast(twi_450m_path)

# Replace -1 values with NA
twi_r[twi_r == -1] <- NA

# Save to NetCDF
terra::writeRaster(
  twi_r,
  twi_450m_mosaic_clean_path,
  filetype  = "GTiff",
  gdal      = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
  overwrite = overwrite,
  datatype  = "FLT4S",
  NAflag    = -9999
)
if (file.exists(twi_450m_mosaic_clean_path)) message("✅ Saved: ", twi_450m_mosaic_clean_path)

