#' Plot Vegetation Height (VEGH) with fixed color scale
#'
#' @param input A SpatRaster object or raster file path
#' @param extent Optional terra::ext() extent.
#'   If provided and smaller than the input raster extent,
#'   the raster will be cropped to their intersecting area.
#' @param title_text Title text for the plot
#' @param text_size Font size
#' @param x_step Number of x-axis breaks
#' @param y_step Number of y-axis breaks
#' @param limits Numeric vector of length 2 specifying color scale limits (optional)
#' @param color_palette Color palette to use (default "batlow")
#' @param color_direction Direction of color palette (1 = normal, -1 = reversed, default -1)
#' @param fill Label for the colorbar (default "m")
#'
#' @return A ggplot2 object
#' @export
plot_vegh <- function(input, extent = NULL, title_text = "Vegetation Height (m)",
                      text_size = 12, x_step = 10, y_step = 10,
                      limits = NULL,
                      color_palette = "batlow",
                      color_direction = -1,
                      fill = "m") {

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

  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = input, maxcell = Inf) +
    scico::scale_fill_scico(
      palette = color_palette,
      direction = color_direction,
      limits = limits,
      na.value = NA,
      oob = scales::squish  # Squish out-of-bounds values to limits
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
