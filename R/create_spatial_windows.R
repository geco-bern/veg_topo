#' Create spatial windows from raster data or data frame
#'
#' Divides spatial data into regular windows and calculates window centroids.
#'
#' @param input Input raster object or data frame
#' @param coord_vars Names of coordinate variables (default: c("lon","lat"))
#' @param value_vars Names of value variables (default: c("twi","vegh"))
#' @param dwin Window size (default: 0.05)
#' @return A data frame with window midpoints and value variables
create_spatial_windows <- function(input,
                                   coord_vars = c("lon", "lat"),
                                   value_vars = c("twi", "vegh"),
                                   dwin = 0.05) {

  # Convert raster to dataframe if needed
  if (inherits(input, "Raster") || inherits(input, "SpatRaster")) {
    suppressWarnings({
      df <- as.data.frame(input, xy = TRUE, na.rm = FALSE)
    })
  } else if (is.data.frame(input)) {
    df <- input
  } else {
    stop("Input must be a raster object or a data frame")
  }

  colnames(df) <- c(coord_vars, value_vars)

  # Create window boundaries
  lon_breaks <- seq(floor(min(df[[coord_vars[1]]])), ceiling(max(df[[coord_vars[1]]])), by = dwin)
  lat_breaks <- seq(floor(min(df[[coord_vars[2]]])), ceiling(max(df[[coord_vars[2]]])), by = dwin)

  # Assign windows and compute midpoints
  df_win <- df |>
    dplyr::mutate(
      ilon = cut(.data[[coord_vars[1]]], breaks = lon_breaks),
      ilat = cut(.data[[coord_vars[2]]], breaks = lat_breaks)
    ) |>
    dplyr::mutate(
      lon_lower = as.numeric(sub("\\((.+),.*", "\\1", ilon)),
      lon_upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", ilon)),
      lat_lower = as.numeric(sub("\\((.+),.*", "\\1", ilat)),
      lat_upper = as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", ilat))
    ) |>
    dplyr::mutate(
      lon_mid = (lon_lower + lon_upper)/2,
      lat_mid = (lat_lower + lat_upper)/2
    ) |>
    dplyr::select(-ilon, -ilat, -lon_lower, -lon_upper, -lat_lower, -lat_upper)

  return(df_win)
}
