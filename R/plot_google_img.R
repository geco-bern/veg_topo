plot_google_img <- function(extent = NULL, title_text = "Google Satellite Map",
                            text_size = 12, x_step = 0.5, y_step = 0.5) {

  # ---- Register Google API key ----
  api_key <- Sys.getenv("GOOGLE_API_KEY")
  register_google(key = api_key)

  # ---- Extract extent boundaries ----
  xmin <- terra::xmin(extent)
  xmax <- terra::xmax(extent)
  ymin <- terra::ymin(extent)
  ymax <- terra::ymax(extent)

  # ---- Compute zoom ----
  earth_circumference <- 40075017
  region_width <- xmax - xmin
  region_width_meters <- region_width * 111000
  zoom_level <- ceiling(log2(earth_circumference / region_width_meters))
  zoom_level <- max(1, min(zoom_level, 21))

  # ---- Fetch map ----
  bbox <- c(left = xmin, bottom = ymin, right = xmax, top = ymax)
  satellite_map <- get_map(location = bbox, source = "google", maptype = "satellite", zoom = zoom_level)

  # ---- Crop to exact extent ----
  # 格式化经度标签
  lon_labels <- function(x) {
    ifelse(x < 0, paste0(abs(x), "°W"),
           ifelse(x > 0, paste0(x, "°E"), "0°"))
  }

  # 格式化纬度标签
  lat_labels <- function(y) {
    ifelse(y < 0, paste0(abs(y), "°S"),
           ifelse(y > 0, paste0(y, "°N"), "0°"))
  }

  p <- ggmap(satellite_map) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
    labs(title = title_text, x = NULL, y = NULL) +
    scale_x_continuous(
      breaks = seq(from = xmin, to = xmax, by = x_step),
      labels = lon_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq(from = ymin, to = ymax, by = y_step),
      labels = lat_labels,
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = text_size),
      axis.text = element_text(size = text_size),
      plot.title = element_text(size = text_size, face = "plain",
                                margin = margin(b = 3)),
      axis.text.x = element_text(
        size = text_size,
        hjust = 0,
        vjust = 1,
        margin = margin(t = 2, b = 2)
      ),
      axis.text.y = element_text(
        size = text_size,
        hjust = 0.5,
        vjust = 0.5,
        angle = 90,
        margin = margin(r = 0, l = 0)
      )
    )

  return(p)
}
