# Merge 3x3 neighboring DEM tiles
# -------helper----------
source(here::here("R/convert_lat.R"))
source(here::here("R/convert_lon.R"))

get_dem_path <- function(lat, lon, file_dir = "") {
  file.path(COP30_dir, sprintf("Copernicus_DSM_10_%s_00_%s_00_DEM.tif",
                               convert_lat(lat), convert_lon(lon)))
}

# ------- Main ----------
merge_dem_neighbors <- function(lat, lon, file_dir = "") {


  dem_file <- get_dem_path(lat, lon, file_dir)

  if(!file.exists(dem_file)) return(NULL)
  dem <- rast(dem_file)

  n_rows <- nrow(dem)
  n_cols <- ncol(dem)

  # Latitudes where neighbor tiles have different longitude resolution!!!!!!
  no_south <- c(-50, -60, -70, -75, -80, -85, 50, 60, 70, 75, 80, 85)
  no_north <- c(-51, -61, -71, -76, -81, -86, 49, 59, 69, 74, 79, 84)

  get_neighbor <- function(dlat, dlon, rows, cols, nrow, ncol) {
    if((dlat == -1 && lat %in% no_south) || (dlat == 1 && lat %in% no_north)) {
      return(matrix(NA_real_, nrow, ncol))
    }
    dem_file <- get_dem_path(lat + dlat, lon + dlon, file_dir)
    if(file.exists(dem_file)) {
      subset <- rast(dem_file)[rows, cols]
      subset <- as.matrix(subset)
      return(matrix(as.numeric(subset), nrow, ncol, byrow=TRUE))
    } else {
      return(matrix(NA_real_, nrow, ncol))
    }
  }

  # Center DEM
  c_matrix <- as.matrix(dem, wide=TRUE)

  # Read all 8 neighbors
  nw <- get_neighbor(1, -1, (n_rows-2):n_rows, (n_cols-2):n_cols, 3, 3)
  n  <- get_neighbor(1, 0,  (n_rows-2):n_rows, 1:n_cols, 3, n_cols)
  ne <- get_neighbor(1, 1,  (n_rows-2):n_rows, 1:3, 3, 3)
  w  <- get_neighbor(0, -1, 1:n_rows, (n_cols-2):n_cols, n_rows, 3)
  e  <- get_neighbor(0, 1,  1:n_rows, 1:3, n_rows, 3)
  sw <- get_neighbor(-1, -1, 1:3, (n_cols-2):n_cols, 3, 3)
  s  <- get_neighbor(-1, 0,  1:3, 1:n_cols, 3, n_cols)
  se <- get_neighbor(-1, 1,  1:3, 1:3, 3, 3)

  # Combine into 3x3 block
  full_matrix <- rbind(
    cbind(nw, n, ne),
    cbind(w,  c_matrix, e),
    cbind(sw, s, se)
  )

  full_raster <- rast(full_matrix)
  ext(full_raster) <- c(
    ext(dem)$xmin - 3 * xres(dem),
    ext(dem)$xmax + 3 * xres(dem),
    ext(dem)$ymin - 3 * yres(dem),
    ext(dem)$ymax + 3 * yres(dem)
  )
  crs(full_raster) <- crs(dem)

  return(full_raster)
}
