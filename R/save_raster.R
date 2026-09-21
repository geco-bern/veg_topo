save_raster <- function(raster, output_file, varname = NULL, overwrite = TRUE) {
  if (!is.null(output_file)) {
    message("Saving...")
    ext <- tolower(tools::file_ext(output_file))

    if (ext == "nc") {
      if (is.null(varname)) varname <- "variable"
      terra::writeCDF(raster, output_file, overwrite = overwrite, varname = varname)
    } else if (ext %in% c("tif", "tiff")) {
      terra::writeRaster(
        raster,
        output_file,
        filetype  = "GTiff",
        gdal      = c("COMPRESS=LZW", "BIGTIFF=YES", "TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256"),
        overwrite = overwrite,
        datatype  = "FLT4S",
        NAflag    = -9999
      )
    } else {
      stop("Only .nc or .tif files supported")
    }

    if (file.exists(output_file)) message("✅ Saved: ", output_file)
  }

  invisible(raster)
}
