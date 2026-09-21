mosaic_tiles_gdal <- function(input_dir, output_file, pattern = "*.tif", overwrite = TRUE, extent = NULL) {

  # Get file list
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found in: ", input_dir)

  message("Found ", length(files), " tiles")

  # Create output directory
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

  # Create file list for gdalbuildvrt
  list_file <- tempfile(fileext = ".txt")
  writeLines(files, list_file)
  on.exit(unlink(list_file))

  # Create VRT
  vrt_file <- tempfile(fileext = ".vrt")
  on.exit(unlink(vrt_file), add = TRUE)

  system(paste("gdalbuildvrt -overwrite -input_file_list",
               shQuote(list_file), shQuote(vrt_file)))

  if (!file.exists(vrt_file)) stop("gdalbuildvrt failed")

  # Build gdal_translate command
  cmd <- "gdal_translate -of COG -co COMPRESS=DEFLATE"

  # Add extent clipping if provided
  if (!is.null(extent)) {
    # Convert SpatExtent to numeric vector using terra::values
    if (inherits(extent, "SpatExtent")) {
      extent <- c(extent$xmin, extent$xmax, extent$ymin, extent$ymax)
    }
    cmd <- paste(cmd, "-projwin", extent[1], extent[4], extent[2], extent[3])
    message("Clipping to: xmin=", extent[1], " xmax=", extent[2],
            " ymin=", extent[3], " ymax=", extent[4])
  }

  cmd <- paste(cmd, shQuote(vrt_file), shQuote(output_file))
  system(cmd)

  # Return raster
  terra::rast(output_file)
}
