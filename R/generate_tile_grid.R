#' Generate Global Tiles Grid
#'
#' Creates a grid of global tiles with specified longitude and latitude steps.
#'
#' @param lon_step Numeric. Longitudinal step size in degrees (default: 30).
#' @param lat_step Numeric. Latitudinal step size in degrees (default: 30).
#' @return A data.frame containing tile names and their geographic extents.

source(here::here("R/make_lon_label.R")) # Helper functions for naming
generate_tile_grid <- function(lon_step = 30, lat_step = 30) {
  # Generate global tile boundaries
  lon_edges <- seq(-180, 180, by = lon_step)
  lat_edges <- seq(-60, 90, by = lat_step)

  # Initialize output data frame
  tile_grid <- data.frame(
    tile_id = character(),
    xmin = numeric(),
    xmax = numeric(),
    ymin = numeric(),
    ymax = numeric(),
    stringsAsFactors = FALSE
  )

  # Generate all possible tiles
  for (lat_i in 1:(length(lat_edges) - 1)) {
    for (lon_j in 1:(length(lon_edges) - 1)) {
      # Create extent object first
      tile_ext <- terra::ext(
        lon_edges[lon_j],
        lon_edges[lon_j + 1],
        lat_edges[lat_i],
        lat_edges[lat_i + 1]
      )

      tile_id <- paste0(
        make_lat_label(lat_edges[lat_i]), "_",
        make_lon_label(lon_edges[lon_j])
      )

      # Use terra's extent accessors
      tile_grid <- rbind(tile_grid, data.frame(
        tile_id = tile_id,
        xmin = terra::xmin(tile_ext),
        xmax = terra::xmax(tile_ext),
        ymin = terra::ymin(tile_ext),
        ymax = terra::ymax(tile_ext),
        stringsAsFactors = FALSE
      ))
    }
  }

  return(tile_grid)
}
