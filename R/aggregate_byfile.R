aggregate_byfile <- function(filename_vegheight, raster_target, outdir){

  # Load the two raster files
  r1 <- rast(filename_vegheight)  # obtained from   https://libdrive.ethz.ch/index.php/s/cO8or7iOe5dT2Rt?path=%2F

  # Aggregate r1 to match r2's resolution
  # The factor is determined by the ratio of resolutions
  fact_x <- (xres(raster_target) / xres(r1))
  fact_y <- (yres(raster_target) / yres(r1))

  # Aggregate using mean (can be changed to other functions like max, min, sum)
  r1_aggregated <- aggregate(r1, fact = c(fact_x, fact_y), fun = mean, na.rm = TRUE)

  # # Optional: Align extents using resample
  # r1_resampled <- resample(r1_aggregated, r2, method = "bilinear")  # Options: "bilinear", "near", etc.

  # create output file name and write to file
  outfilnam <- paste0(outdir, str_remove(basename(filename_vegheight), ".tif"), "_15arcsec.nc")
  message(paste("Writing to file", outfilnam, "..."))
  writeCDF(r1_aggregated, outfilnam, overwrite = TRUE)
}
