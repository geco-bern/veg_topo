# Mosaic raster tiles with optional resampling/masking
#
# Args:
#   input_dir    Directory with input tiles
#   output_file  Output path (NULL = do not save)
#   target_grid  Template raster or path (for resample/mask)
#   if_resample  Apply resampling to target grid
#   if_mask      Apply masking using target grid
#   pattern      File pattern for tiles (e.g. "*.nc")
#   overwrite    Overwrite existing file
#   crs          CRS to assign (optional)
#
# Returns:
#   SpatRaster (invisible)
source(here::here("R/save_raster.R"))
mosaic_tiles <- function(input_dir,
                         output_file = NULL,
                         target_grid = NULL,
                         if_resample = FALSE,
                         if_mask = FALSE,
                         pattern = "*.nc",
                         overwrite = TRUE,
                         crs = NULL,
                         varname = "band",
                         if_crop = FALSE) {

  # ===========================================================================
  # 1) Input check & file listing
  # ===========================================================================
  tile_paths <- fs::dir_ls(path = input_dir, glob = pattern)

  if (length(tile_paths) == 0) {
    stop("No matching files found in: ", input_dir)
  }

  # ===========================================================================
  # 2) Target grid (optional)
  # ===========================================================================
  if (!is.null(target_grid)) {
    r_tar <- if (is.character(target_grid)) terra::rast(target_grid)[[1]] else target_grid[[1]]
  }

  # ===========================================================================
  # 3) Read tiles & mosaic
  # ===========================================================================
  message("Mosaicing...")
  tile_rasters <- lapply(tile_paths, terra::rast)
  r_out <- do.call(terra::merge, unname(tile_rasters))

  # Assign CRS if provided
  if (!is.null(crs)) {
    terra::crs(r_out) <- crs
  }

  # ===========================================================================
  # 4) Resample / mask (optional)
  # ===========================================================================
  if (if_resample && !is.null(target_grid)) {
    message("Resampling...")
    r_out <- terra::resample(r_out, r_tar, method = "bilinear")
  }

  if (if_mask && !is.null(target_grid)) {
    message("Masking...")
    r_out <- terra::mask(r_out, r_tar)
    }

  if (if_crop && !is.null(target_grid)) {
    message("Cropping...")
    r_out <- terra::crop(r_out, r_tar)
  }
  # ===========================================================================
  # 5) Write output (optional)
  # ===========================================================================
  save_raster(r_out, output_file, varname)

  # ===========================================================================
  # 6) Clean up & return
  # ===========================================================================
  rm(tile_rasters); gc()

  return(invisible(r_out))
}
