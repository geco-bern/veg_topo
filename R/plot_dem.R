
plot_dem <- function(input, extent = NULL, title_text = "Elevation (km)", fill = "km",
                     limits = NULL, text_size = 12, x_step = 10, y_step = 10) {

  # ---- Load raster ----
  if (is.character(input)) input <- terra::rast(input)
  if (!inherits(input, "SpatRaster")) stop("Input must be a SpatRaster or valid file path.")

  # ---- Handle extent and optional cropping ----
  if (is.null(extent)) {
    extent <- terra::ext(input)
  } else if (!inherits(extent, "SpatExtent")) {
    stop("`extent` must be a SpatExtent object from terra::ext().")
  } else {
    # Crop raster if extent is smaller
    area_in <- (terra::xmax(input) - terra::xmin(input)) * (terra::ymax(input) - terra::ymin(input))
    area_ex <- (terra::xmax(extent) - terra::xmin(extent)) * (terra::ymax(extent) - terra::ymin(extent))

    if (area_ex < area_in) {
      cropped <- terra::crop(input, extent)
      if (all(is.na(terra::values(cropped)))) stop("Extent does not intersect raster.")
      input <- cropped
      message("Raster cropped to specified extent.")
    }
  }

  # ---- Extract extent boundaries ----
  xmin <- terra::xmin(extent)
  xmax <- terra::xmax(extent)
  ymin <- terra::ymin(extent)
  ymax <- terra::ymax(extent)

  # Compute value range for color scale (if limits not provided)
  if (is.null(limits)) {
    vmin <- terra::global(input, "min", na.rm = TRUE)[1, 1] |> as.numeric()
    vmax <- terra::global(input, "max", na.rm = TRUE)[1, 1] |> as.numeric()
    limits <- c(vmin, vmax)
  }

  # ---- Plot ----
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = input, maxcell = Inf) +
    scale_fill_gradientn(
      colors = terrain.colors(255),
      # direction = 1,
      na.value = NA,
      limits = limits,
      labels = function(x) x / 1000
    ) +
    ggplot2::labs(
      title = title_text,
      fill = fill,
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(from = xmin, to = xmax, by = x_step),
      expand = c(0, 0),
      limits = c(xmin, xmax)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(from = ymin, to = ymax, by = y_step),
      expand = c(0, 0),
      limits = c(ymin, ymax)
    ) +
    ggplot2::theme_bw(base_size = text_size)

  return(p)
}
