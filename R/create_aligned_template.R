source(here::here("R/get_lonlat_extent.R"))
create_aligned_template <- function(input, res_out = 0.05, crs_out = "EPSG:4326", snap_step = NULL, trim_input = FALSE) {
  # ---- Handle different input types ----
  if (inherits(input, "SpatExtent")) {
    # Input is already an extent
    e <- input
    # Need CRS info - if not provided, use crs_out
    if (is.null(crs_out)) {
      stop("When input is an extent, crs_out must be specified")
    }
    input_crs <- crs_out
  } else {
    # Original logic for raster inputs
    if (is.character(input)) input <- terra::rast(input)
    if (!inherits(input, "SpatRaster")) stop("Input must be a SpatRaster, valid file path, or SpatExtent.")

    if (trim_input ) input <- trim(input)

    # Get input raster extent
    e <- if(terra::crs(input) != "EPSG:4326") get_lonlat_extent(input) else terra::ext(input)
    input_crs <- terra::crs(input)
  }

  # ---- Helper function to determine decimal places ----
  get_decimal_places <- function(x) {
    if (x == 0) return(0)
    x_str <- format(x, scientific = FALSE)
    if (!grepl("\\.", x_str)) return(0)
    nchar(strsplit(x_str, "\\.")[[1]][2])
  }

  # ---- Set alignment step and rounding precision ----
  if(is.null(snap_step)) snap_step <- res_out

  # ---- Avoid floating-point errors using integer arithmetic ----
  scale_factor <- 1 / snap_step  # 0.00025 -> 4000

  xmin_aligned <- floor(round(e$xmin * scale_factor, 0)) / scale_factor
  xmax_aligned <- ceiling(round(e$xmax * scale_factor, 0)) / scale_factor
  ymin_aligned <- floor(round(e$ymin * scale_factor, 0)) / scale_factor
  ymax_aligned <- ceiling(round(e$ymax * scale_factor, 0)) / scale_factor

  # Compute number of rows and columns to cover the aligned extent
  ncols <- round((xmax_aligned - xmin_aligned) / res_out)
  nrows <- round((ymax_aligned - ymin_aligned) / res_out)

  # Use input CRS if crs_out not provided and input is raster
  if (is.null(crs_out)) {
    if (exists("input_crs") && !inherits(input, "SpatExtent")) {
      crs_out <- input_crs
    } else {
      crs_out <- "EPSG:4326"  # Default fallback
    }
  }

  # Create the template raster
  template <- terra::rast(
    nrows = nrows,
    ncols = ncols,
    xmin = xmin_aligned,
    xmax = xmax_aligned,
    ymin = ymin_aligned,
    ymax = ymax_aligned,
    crs = crs_out
  )

  return(template)
}
