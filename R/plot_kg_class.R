#' Plot Köppen–Geiger Climate Classification Raster
#'
#' This function plots a Köppen–Geiger climate classification raster using ggplot2 and tidyterra.
#' It supports optional cropping by extent and automatically filters the legend to only show
#' climate classes present in the raster.
#'
#' @param input SpatRaster or file path to raster
#' @param legend_file Path to a legend text file defining class codes, names, and colors
#' @param extent Optional SpatExtent to crop the raster
#' @param title_text Plot title (default: "Köppen–Geiger Class")
#' @param text_size Base text size for plot
#' @param x_step Interval for x-axis ticks
#' @param y_step Interval for y-axis ticks
#'
#' @return ggplot object
plot_kg_class <- function(input, legend_file = NULL, extent = NULL, title_text = "Köppen–Geiger Class",
                          text_size = 12, x_step = 10, y_step = 10) {

  # ---- Load raster ----
  if (is.character(input)) input <- terra::rast(input)
  if (!inherits(input, "SpatRaster")) stop("Input must be a SpatRaster or valid file path.")

  # ---- Handle extent and crop raster if needed ----
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

  # ---- Load legend ----
  legend_lines <- readLines(legend_file)
  class_lines <- legend_lines[grep("^\\s+\\d+:", legend_lines)]

  kg_simple <- data.frame(
    code = as.numeric(sub(":.*", "", gsub(" ", "", class_lines))),
    class = sapply(strsplit(class_lines, "\\s+"), function(x) x[3]),
    color_hex = sapply(regmatches(class_lines, regexpr("\\[[0-9 ]+\\]", class_lines)),
                       function(x) {
                         rgb_vec <- as.numeric(strsplit(gsub("\\[|\\]", "", x), " ")[[1]])
                         rgb(rgb_vec[1], rgb_vec[2], rgb_vec[3], maxColorValue = 255)
                       }),
    stringsAsFactors = FALSE
  )

  # ---- Filter legend to actual raster values ----
  actual_values <- na.omit(unique(terra::values(input)))
  actual_values <- as.character(actual_values)

  kg_simple$code_char <- as.character(kg_simple$code)
  existing_legend <- kg_simple[kg_simple$code_char %in% actual_values, ]

  # ---- Convert raster to factor for plotting ----
  kg_factor <- as.factor(input)

  # ---- Plot raster ----
  p <- ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = kg_factor, maxcell = Inf) +
    scale_fill_manual(
      name = NULL,
      values = setNames(existing_legend$color_hex, existing_legend$code_char),
      labels = setNames(existing_legend$class, existing_legend$code_char),
      guide = guide_legend(keywidth = 0.8),
      na.value = NA
    ) +
    ggplot2::labs(
      title = title_text,
      fill = NULL
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
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = text_size * 0.9),
      legend.title = ggplot2::element_text(size = text_size),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, -8),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = text_size * 0.8,
        hjust = 0.5,
        vjust = 0.5,
        margin = margin(t = 2, b = 2),
      ),
      axis.text.y = ggplot2::element_text(
        size = text_size * 0.8,
        hjust = 0.5,
        vjust = 0.5,
        margin = margin(r = 0, l = 2)
      ),
      panel.spacing = unit(0, "cm"),
      panel.border = ggplot2::element_rect(linewidth = 0.5, fill = NA),
      plot.margin = margin(0, 0, 0, 0),
      plot.title = ggplot2::element_text(
        size = text_size * 1.2,
        face = "plain",
        margin = margin(b = 0)
      ),
      plot.title.position = "panel"
    )

  return(p)
}

