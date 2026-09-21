# ~ 1 min on UBELIX
# -------------------- Set Up --------------------------------------------------
library(terra)
library(dplyr)
library(arrow)

# Load custom functions
source(here::here("R/config.R"))
source(here::here("R/create_aligned_template.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/convert_lat.R"))
source(here::here("R/convert_lon.R"))

# Create output directory
if (!dir.exists(rf_sample_data_tiles_dir)) {
  dir.create(rf_sample_data_tiles_dir, recursive = TRUE)
  message("Directory created: ", rf_sample_data_tiles_dir)
}

# -------------------- Prepare Valid Windows ------------------------------------
# Load correlation rasters
r_h_twi <- rast(cor_twi_vegh_mosaic_file)
r_h_rin <- rast(r_H_R_5km_path)
fused <- rast(fused_5km_file)

# Identify valid windows (both correlations available)
valid_win <- !is.na(r_h_twi) & !is.na(r_h_rin) & fused < 0.05

# Statistics
true_count <- global(valid_win, "sum", na.rm = TRUE)[[1]]
message(sprintf("Valid windows (TRUE): %s", format(true_count, big.mark = ",")))
# Valid windows (TRUE): 4,853,773

# Save
terra::writeRaster(valid_win, valid_win_path, overwrite = TRUE)
if(file.exists(valid_win_path)) message(sprintf("Saved: %s", valid_win_path))
# -------------------- Create 1° Tiles -----------------------------------------
# Create 1° template
align_1_deg <- create_aligned_template(valid_win, res_out = 1)

# Aggregate to 1° tiles
# Only consider pixels with cover fraction > 0.1 to eliminate border effects
grid_1_deg <- raster_preprocess_save(
  input        = valid_win,
  output       = NULL,
  target       = align_1_deg,
  varname      = "valid_tile",
  if_zonal     = TRUE,
  fun          = function(values, coverage_fractions) {
    valid_indices <- !is.na(values) & coverage_fractions > 0.1
    if (any(valid_indices)) {
      any(values[valid_indices] == 1, na.rm = TRUE)
    } else {
      FALSE
    }
  },
  if_aggregate = FALSE,
  if_resample  = FALSE,
  if_return_raster = TRUE
)

# Summary
cat("Valid tiles summary:\n")
print(table(values(grid_1_deg), useNA = "no"))
# FALSE  TRUE
# 34634 17206
# Total tiles to process: 17206
# -------------------- Create Tile Dataframe ------------------------------------
rf_tiles <- as.data.frame(grid_1_deg, xy = TRUE) |>
  filter(lyr.1 == TRUE) |>
  transmute(
    xmin = floor(x - 0.5),
    xmax = ceiling(x + 0.5),
    ymin = floor(y - 0.5),
    ymax = ceiling(y + 0.5)
  ) |>
  mutate(
    reg_id = paste0(convert_lat(ymin), "_", convert_lon(xmin))
  )

message(sprintf("Total tiles to process: %d", nrow(rf_tiles)))

# Save
arrow::write_parquet(rf_tiles, rf_tiles_path)
if(file.exists(rf_tiles_path)) message(sprintf("Saved: %s", rf_tiles_path))
