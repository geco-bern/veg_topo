create_conical_dem <- function(center_lon, center_lat, slope_angle_deg,
                               output_file = NULL, window_size = 0.1) {

  # Convert slope angle to elevation (meters)
  slope_rate <- tan(slope_angle_deg * pi / 180)
  res <- 15 / 3600  # degrees
  n <- floor(window_size / res)
  if(n %% 2 == 0) n <- n + 1

  # Create equally spaced coordinates
  half_dist_deg <- window_size / 2
  lon <- seq(center_lon - half_dist_deg, center_lon + half_dist_deg, length.out = n)
  lat <- seq(center_lat - half_dist_deg, center_lat + half_dist_deg, length.out = n)

  # Create grids
  lon_grid <- matrix(lon, n, n, byrow = TRUE)
  lat_grid <- matrix(lat, n, n, byrow = FALSE)

  # Calculate distance from each point to center (degrees)
  dist_deg <- sqrt((lon_grid - center_lon)^2 + (lat_grid - center_lat)^2)

  # Convert to meters (using center latitude)
  meters_per_deg_lat <- 111320
  meters_per_deg_lon <- 111320 * cos(center_lat * pi / 180)

  # Calculate distance (meters)
  dist_x_m <- (lon_grid - center_lon) * meters_per_deg_lon
  dist_y_m <- (lat_grid - center_lat) * meters_per_deg_lat
  dist_meters <- sqrt(dist_x_m^2 + dist_y_m^2)

  # Calculate circle radius (meters)
  radius_m <- (window_size / 2) * meters_per_deg_lon  # Using longitudinal distance as radius

  # Calculate elevation: cone inside circle, 0 outside
  elev <- matrix(0, n, n)
  inside_circle <- dist_meters <= radius_m

  # Elevation inside circle (conical shape)
  max_height_m <- slope_rate * radius_m
  elev[inside_circle] <- max_height_m - slope_rate * dist_meters[inside_circle]
  elev[elev < 0] <- 0

  # Create raster
  r <- terra::rast(nrows = n, ncols = n,
                   xmin = min(lon), xmax = max(lon),
                   ymin = min(lat), ymax = max(lat),
                   crs = "EPSG:4326")
  terra::values(r) <- as.vector(elev)

  if(!is.null(output_file)) terra::writeRaster(r, output_file, overwrite = TRUE)
  return(r)
}
