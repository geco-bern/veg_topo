#' Plot Correlation between VEGH and R (Radiation)
#'
#' @param input A SpatRaster object or a raster file path
#' @param extent Optional terra::ext() extent.
#'   If provided and smaller than the input raster extent,
#'   the raster will be cropped to their intersecting area.
#' @param text_size Font size
#' @param x_step Number of x-axis step
#' @param y_step Number of y-axis step
#' @return A ggplot2 object
#' @export
plot_r_H_R <- function(input, extent = NULL, title_text = "Pearson's r (H～Rᵢₙ)",fill = "",
                       text_size = 12, x_step = 10, y_step = 10) {
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

  # ---- Plot ----
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = input, maxcell = Inf) +
    scico::scale_fill_scico(
      palette = "vik",
      direction = -1,
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.5),
      midpoint = 0,
      # name = expression(r[H*","*R]),
      na.value = NA
    ) +
    ggplot2::labs(
      title = title_text,
      fill = fill,
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(from = xmin, to = xmax, by = x_step),
      # limits = c(xmin, xmax),
      expand = expansion(mult = 0.0001)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(from = ymin, to = ymax, by = y_step),
      # limits = c(ymin, ymax),
      expand = expansion(mult = 0.0001)
    ) +
    ggplot2::theme_bw(base_size = text_size)

  return(p)
}
