
# ~ 120 min on UBELIX
# ------------- Setup-----------------------------------------------------------

library(terra)
library(fs)

source(here::here("R/config.R"))
source(here::here("R/mosaic_tiles_gdal.R"))

ext <- ext(rast(twi_450m_mosaic_clean_path))
message("Output raster will be cliped to", ext)

# ------------- Mosaicing DEM------------------------------------------------------
tictoc::tic()

message("Mosaicing DEM...")
dem_mosaic <- mosaic_tiles_gdal(
  input_dir   = dem_450m_tiles_dir,
  output_file = dem_450m_mosaic_path,
  extent = ext)

tictoc::toc()

# ------------- Mosaicing Slope------------------------------------------------------
tictoc::tic()

message("Mosaicing Slope...")
slope_mosaic <- mosaic_tiles_gdal(
  input_dir   = slope_450m_tiles_dir,
  output_file = slope_450m_mosaic_path,
  extent = ext)

tictoc::toc()

# ------------- Mosaicing Aspect------------------------------------------------------
tictoc::tic()

message("Mosaicing Aspect...")
aspect_mosaic <- mosaic_tiles_gdal(
  input_dir   = aspect_450m_tiles_dir,
  output_file = aspect_450m_mosaic_path,
  extent = ext)

tictoc::toc()
# ------------- Mosaicing Rin------------------------------------------------------

tictoc::tic()
message("Mosaicing Rin...")
rin_mosaic <- mosaic_tiles_gdal(
  input_dir   = rin_450m_tiles_dir,
  output_file = rin_450m_mosaic_path,
  extent = ext)
tictoc::toc()
