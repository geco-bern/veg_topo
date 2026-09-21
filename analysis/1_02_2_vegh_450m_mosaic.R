# ～ 15 min： on UBELIX
# rely on the output from ~/veg_topo/analysis/1_01_1_twi_450m_clean.R !!!!

# --------------- Setup --------------------------------------------------------
library(terra)

source(here::here("R/config.R"))
source(here::here("R/mosaic_tiles.R"))

# --------------- Mosaic, resample and save ------------------------------------
tictoc::tic()

message("Vegetation height...")
vegh_450m_mosaic <- mosaic_tiles(
  input_dir   = vegh_450m_tiles_dir,
  output_file = vegh_450m_mosaic_path, # Save as GeoTiff
  pattern = "*_Map_to450m.nc",
  target_grid = twi_450m_mosaic_clean_path,
  if_resample = TRUE
  )

message("Fraction of vegetated area...")
fveg_450m_mosaic <- mosaic_tiles(
  input_dir   = vegh_450m_tiles_dir,
  output_file = fveg_real_450m_mosaic_path,  # Save as GeoTiff
  pattern = "*_Map_to450m_fveg.nc",
  target_grid = twi_450m_mosaic_clean_path,
  if_resample = TRUE
  )

tictoc::toc()

# ----------------- Check result (Optional) -----------------
# r <- rast(vegh_450m_mosaic_path)
# plot(r)

# ----------------- Delete intermedia files  (Optional) -----------------
# unlink(vegh_450m_tiles_dir, recursive = TRUE)
