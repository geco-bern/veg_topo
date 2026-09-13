source(here::here("R/save_raster.R"))
df_to_raster <- function(df_calc,
                         x = "lon",
                         y = "lat",
                         value = "rin",
                         template_raster = NULL,
                         output_file = NULL,
                         varname = "band",
                         overwrite = TRUE,
                         return_raster = TRUE) {

  # Build raster from xyz
  r_out <- terra::rast(df_calc[, c(x, y, value)],
                       type = "xyz",
                       crs = crs(template_raster))

  # Extend and mask
  r_out <- terra::extend(r_out, template_raster)

  # Save if output_file is provided

  if (!is.null(output_file)) save_raster(r_out, output_file)

  # Return based on option
  if (return_raster) return(invisible(r_out))
  return(invisible(output_file))
}
