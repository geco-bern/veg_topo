#' Get the extent of a raster in WGS84 longitude/latitude coordinates
#'
#' Reprojects a dense grid of points from the raster's CRS to WGS84 and returns
#' the bounding box. Using a dense grid (not just corners) captures the true
#' extent, accounting for projection distortion and edge curvature.
#'
#' @param r A SpatRaster object
#' @param n Number of sample points along each axis (default: 50).
#'          Total points = n * n.
#'
#' @return A SpatExtent object in WGS84 (EPSG:4326)
#'
#' @examples
#' r <- rast("some_projected_raster.tif")
#' lonlat_ext <- get_lonlat_extent(r)
get_lonlat_extent <- function(r, n = 50) {
  e <- ext(r)

  # Create a grid of sample points across the raster (not just 4 corners)
  # This captures the true shape after reprojection
  x_pts <- seq(e$xmin, e$xmax, length.out = n)
  y_pts <- seq(e$ymin, e$ymax, length.out = n)
  pts <- expand.grid(x = x_pts, y = y_pts)

  # Transform all points from raster projection to WGS84 geographic coordinates
  sf_pts <- sf::st_as_sf(pts, coords = c("x", "y"), crs = crs(r))
  sf_pts_wgs84 <- sf::st_transform(sf_pts, crs = 4326)
  coords <- sf::st_coordinates(sf_pts_wgs84)

  # Return the bounding box in degrees (lon/lat)
  # Using dense sampling ensures we don't miss the true min/max due to projection distortion
  return(ext(min(coords[,1]), max(coords[,1]),
             min(coords[,2]), max(coords[,2])))
}
