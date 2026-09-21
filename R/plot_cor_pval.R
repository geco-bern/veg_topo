#' Plot correlation p-values as categorized raster
#'
#' @param input SpatRaster object or path to raster file containing p-values
#' @param extent Optional SpatExtent to crop the raster
#' @param title_text Plot title
#' @param text_size Base font size
#' @param x_step Interval for x-axis ticks
#' @param y_step Interval for y-axis ticks
#' @return ggplot2 object
#' @export

plot_cor_pval <- function(input, extent = NULL,
                          title_text = "Pearson's p-value (H~TWI)",
                          fill_text = "",
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

  # ---- 获取图层名称 ----
  # 获取第一个图层的名称用于美学映射
  layer_name <- names(input)[1]

  breaks_sqrt <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1)
  labels_sqrt <- c("0.01", "0.05", "0.10", "0.25", "0.50", "0.75", "1.00")
  # ---- Plot using tidyterra ----
  p <- ggplot() +
    tidyterra::geom_spatraster(
      data = input,
      maxcell = Inf
    ) +
    scico::scale_fill_scico(
      palette = "batlowK",
      direction = -1,
      na.value = NA,
      trans = "sqrt",  # 平方根变换
      breaks = breaks_sqrt,
      labels = labels_sqrt,
      guide = guide_colorbar(
        barwidth = 0.8,
        barheight = 6,
        title = NULL,
        direction = "vertical"
      )
    ) +
    ggplot2::labs(
      title = title_text,
      fill = fill_text,
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(from = xmin, to = xmax, by = x_step),
      limits = c(xmin, xmax),
      expand = expansion(mult = 0.0001)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(from = ymin, to = ymax, by = y_step),
      limits = c(ymin, ymax),
      expand = expansion(mult = 0.0001)
    ) +
    ggplot2::theme_bw(base_size = text_size)

  return(p)
}
