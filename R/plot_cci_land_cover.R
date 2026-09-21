plot_cci_land_cover <- function(input, extent = NULL,
                                title_text = "Land Cover Type",
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

  # Load raster layer
  lccs_class <- input[["lccs_class"]]

  # Crop raster to extent
  landcover_crop <- terra::crop(lccs_class, extent)

  # 检查裁剪后的数据是否有效
  if (terra::ncell(landcover_crop) == 0) {
    stop("Cropped raster has no cells.")
  }

  # Convert to categorical raster and assign labels/colors
  landcover_crop <- terra::as.factor(landcover_crop)

  # 创建颜色映射数据框
  color_df <- data.frame(
    value = c(0, 10, 11, 12, 20, 30, 40, 50, 60, 61, 62,
              70, 71, 72, 80, 81, 82, 90, 100, 110, 120,
              121, 122, 130, 140, 150, 151, 152, 153, 160,
              170, 180, 190, 200, 201, 202, 210, 220),
    label = c("No Data", "Cropland, rainfed", "Herbaceous cover", "Tree or shrub cover",
              "Cropland, irrigated or post-flooding",
              "Mosaic cropland >50% / natural vegetation <50%",
              "Mosaic natural vegetation >50% / cropland <50%",
              "Tree broadleaf evergreen", "Tree broadleaf deciduous",
              "Tree broadleaf deciduous closed", "Tree broadleaf deciduous open",
              "Tree needleleaf evergreen", "Tree needleleaf evergreen closed",
              "Tree needleleaf evergreen open", "Tree needleleaf deciduous",
              "Tree needleleaf deciduous closed", "Tree needleleaf deciduous open",
              "Tree mixed leaf type", "Mosaic tree/shrub >50%",
              "Mosaic herbaceous >50%", "Shrubland",
              "Evergreen shrubland", "Deciduous shrubland", "Grassland",
              "Lichens and mosses", "Sparse vegetation <15%", "Sparse tree <15%",
              "Sparse shrub <15%", "Sparse herbaceous <15%",
              "Tree flooded fresh/brackish", "Tree flooded saline",
              "Shrub/herb flooded", "Urban areas", "Bare areas",
              "Consolidated bare", "Unconsolidated bare",
              "Water bodies", "Snow and ice"),
    color = c("#000000", "#FFFF64", "#FFFF64", "#FFFF00", "#AAF0F0", "#DCF064", "#C8C864",
              "#006400", "#00A000", "#00A000", "#AAC800", "#003C00", "#003C00", "#005000",
              "#285000", "#285000", "#326400", "#788000", "#8CA000", "#BE9600", "#966400",
              "#966400", "#966400", "#FFB432", "#FFDCD6", "#FFEBAF", "#FFC864", "#FFD278",
              "#FFEBAF", "#00785A", "#009678", "#00DC82", "#C31400", "#FFF5D7", "#DCDCDC",
              "#FFF5D7", "#0046C8", "#FFFFFF")
  )

  # 创建颜色和标签映射
  color_map <- setNames(color_df$color, as.character(color_df$value))
  label_map <- setNames(as.character(color_df$value), as.character(color_df$value))

  # 只保留裁剪后实际存在的类别
  unique_vals <- unique(terra::values(landcover_crop))
  unique_vals <- unique_vals[!is.na(unique_vals)]

  if (length(unique_vals) > 0) {
    # 过滤颜色映射，只包含实际存在的类别
    color_map <- color_map[names(color_map) %in% as.character(unique_vals)]
    label_map <- label_map[names(label_map) %in% as.character(unique_vals)]
  }

  # Plot
  p <- ggplot() +
    tidyterra::geom_spatraster(data = landcover_crop) +
    scale_fill_manual(
      values = color_map,
      labels = label_map,
      name = "",
      guide = guide_legend(keywidth = 0.8, keyheight = 0.4),
      na.translate = FALSE  # 不显示NA值的图例
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Longitude",
      y = "Latitude",
      fill = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(from = xmin, to = xmax, by = x_step),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(from = ymin, to = ymax, by = y_step),
      expand = c(0, 0)
    ) +
    ggplot2::coord_sf(
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      expand = FALSE
    ) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = text_size * 0.9),
      legend.title = ggplot2::element_text(size = text_size),
      axis.title = ggplot2::element_text(size = text_size),
      axis.text = ggplot2::element_text(size = text_size * 0.9),
      plot.title = ggplot2::element_text(size = text_size * 1.2, face = "bold"),
      plot.title.position = "panel"
    )

  return(p)
}
