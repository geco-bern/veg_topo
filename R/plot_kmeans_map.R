plot_kmeans_map <- function(input, extent = NULL, title_text = "K-means cluster map (k=8)",
                            highlight_cluster = NULL,
                            text_size = 12, x_step = 30, y_step = 30) {

  # Reorder cluster_value and cluster_labels with an defined order
  # load(here::here("data/cluster_data.RData")) # cluster_values, cluster_labels

  # fill_color for dry to wet cluster
  fill_colors <- c(
    "#E78AC3",  # Pink - Arid
    "#FC8D62",  # Orange - Semi-arid
    "#FFD92F",  # Yellow - Semi-arid
    "#E5C494",  # Light brown - Dry-sub-humid
    "#B3B3B3",  # Gray - Humid
    "#66C2A5", # Blue-green - Humid
    "#8DA0CB",  # Blue - Humid
    "#A6D854"   # Green - Humid
  )

  # Remove names from colors
  fill_colors <- unname(fill_colors)

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

  # Get actual values present in the (possibly cropped) raster
  raster_vals <- terra::values(input) |> na.omit() |> unique() |> sort()

  # ---- Handle cluster labels ----
  # 确保cluster_labels正确命名（使用cluster_values的顺序）
  names(cluster_labels) <- as.character(cluster_values)

  # Get the actual categories present after cropping
  if (!is.null(cluster_labels)) {
    # Filter cluster_labels to only include those present in the cropped raster
    present_clusters <- as.character(raster_vals)

    # 关键修改：保持cluster_labels的原始顺序
    # cluster_labels是按照cluster_values顺序命名的
    # 所以我们需要保持这个顺序

    # 获取所有存在的聚类（按cluster_values顺序）
    all_clusters <- as.character(cluster_values)
    existing_clusters <- all_clusters[all_clusters %in% present_clusters]

    # 获取对应的标签（保持cluster_values的顺序）
    present_labels <- cluster_labels[existing_clusters]

    if (length(present_labels) > 0) {
      # Create mapping between numeric values and labels
      # 保持cluster_values中的顺序，而不是数字大小顺序
      value_to_label <- setNames(
        as.character(present_labels),
        as.numeric(names(present_labels))
      )

      # Get final levels in correct order（按照cluster_values的顺序）
      final_levels <- unname(value_to_label)
    } else {
      # If no cluster labels match, use numeric values
      value_to_label <- setNames(as.character(raster_vals), raster_vals)
      final_levels <- as.character(raster_vals)
    }
  } else {
    # Use raster unique values directly if no labels provided
    value_to_label <- setNames(as.character(raster_vals), raster_vals)
    final_levels <- as.character(raster_vals)
  }

  # ---- Create color mappings ----
  # Match colors to present levels
  if (!is.null(cluster_labels) && exists("present_labels")) {
    # 获取当前存在的聚类数值
    present_nums <- as.numeric(names(present_labels))

    # 确保按照cluster_values的顺序分配颜色
    # cluster_values: 2 3 1 5 4 6 8 7
    # fill_colors: 1st, 2nd, 3rd, 4th, 5th, 6th, 7th, 8th
    # 所以：聚类2 -> 颜色1, 聚类3 -> 颜色2, 聚类1 -> 颜色3, 聚类5 -> 颜色4...

    # 查找每个聚类值在cluster_values中的位置
    color_positions <- match(present_nums, cluster_values)
    present_colors <- fill_colors[color_positions]

  } else {
    # For numeric values, use colors based on sorted values
    color_index <- match(raster_vals, 1:length(fill_colors))
    present_colors <- fill_colors[color_index]
  }

  # Create named vectors for scale_* functions
  # Ensure colors are named by the labels (final_levels)
  names(present_colors) <- final_levels

  # ---- Prepare data for plotting ----
  # Convert raster to dataframe for better control
  raster_df <- as.data.frame(input, xy = TRUE, na.rm = TRUE)
  colnames(raster_df) <- c("x", "y", "value")

  # Debug: print unique values in dataframe
  message("Unique values in raster_df: ", paste(unique(raster_df$value), collapse = ", "))
  message("Value to label mapping: ", paste(names(value_to_label), "->", value_to_label, collapse = ", "))

  # Convert value to factor with proper labels
  # IMPORTANT: Ensure we use the same factor levels as in value_to_label
  raster_df$value_factor <- factor(
    raster_df$value,
    levels = as.numeric(names(value_to_label)),
    labels = unname(value_to_label)
  )

  # Debug: print factor levels
  message("Factor levels: ", paste(levels(raster_df$value_factor), collapse = ", "))
  message("Factor values in data: ", paste(head(levels(raster_df$value_factor)[raster_df$value_factor]), collapse = ", "))

  # ---- Handle highlighting ----
  if (!is.null(highlight_cluster)) {
    # Find the corresponding label for the highlighted cluster
    # First try to get from cluster_labels
    highlight_label <- NULL

    if (is.numeric(highlight_cluster) || (is.character(highlight_cluster) && grepl("^[0-9]+$", highlight_cluster))) {
      # If it's a number (or numeric string), look it up in cluster_labels
      highlight_num <- as.numeric(highlight_cluster)
      if (highlight_num %in% as.numeric(names(cluster_labels))) {
        highlight_label <- as.character(cluster_labels[as.character(highlight_num)])
      }
    } else {
      # If it's already a label, check if it exists in final_levels
      if (highlight_cluster %in% final_levels) {
        highlight_label <- highlight_cluster
      }
    }

    # Debug highlight
    message("Highlight cluster input: ", highlight_cluster)
    message("Highlight label found: ", highlight_label)

    if (!is.null(highlight_label)) {
      # Create alpha column based on highlight
      raster_df$alpha <- ifelse(as.character(raster_df$value_factor) == highlight_label, 1, 0)
      message("Alpha values: ", paste(table(raster_df$alpha), collapse = ", "))
    } else {
      warning("Highlight cluster not found in data. No highlighting applied.")
      raster_df$alpha <- 1
    }
  } else {
    raster_df$alpha <- 1
  }

  # ---- Create ggplot ----
  p <- ggplot() +
    geom_tile(
      data = raster_df,
      aes(x = x, y = y, fill = value_factor, alpha = alpha)
    ) +
    scale_fill_manual(
      values = present_colors,
      name = "Group",
      na.value = NA,
      na.translate = FALSE,
      guide = guide_legend(
        keywidth = unit(0.7, "in"),
        keyheight = unit(0.1, "in"),
        # override.aes = list(size = 5),
        title.position = "left",
        label.position = "bottom",
        nrow = 1
      )
    ) +
    scale_alpha_identity(guide = "none") +  # Use identity scale for alpha
    labs(title = title_text) +
    scale_x_continuous(
      breaks = seq(from = xmin, to = xmax, by = x_step),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(from = ymin, to = ymax, by = y_step),
      expand = c(0, 0)
    ) +
    coord_sf(
      xlim = c(xmin, xmax),
      ylim = c(ymin, ymax),
      expand = FALSE,
      clip = "off"
    ) +
    theme_bw(base_size = text_size) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = text_size * 0.9),
      legend.title = element_text(size = text_size),
      legend.background = element_blank(),  # 图例整体背景透明
      legend.box.background = element_blank(),  # 图例框背景透明
      axis.title = element_text(size = text_size),
      axis.text = element_text(size = text_size * 0.9),
      plot.title = element_text(size = text_size * 1.2, face = "bold"),
      plot.title.position = "panel"
    )

  return(p)
}
