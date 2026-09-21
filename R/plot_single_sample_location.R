
#-------------------------------------------------------------
# Function: plot_single_sample_location
# Description:
#   Given a coordinate (lon, lat), this function selects the
#   nearest predefined 90°×90° window center and plots the
#   corresponding region on a world map.
#-------------------------------------------------------------

# library(sf)
# library(ggplot2)
# library(rnaturalearth)

plot_single_sample_location <- function(lon, lat,  title_text = "Location", text_size = 12, point_color = "darkred", point_size = 1) {

  # Define custom window centers
  centers_df <- data.frame(
    lon_c = c(-135, -105, -75,   -0,  30,   60,  90,   135,   -60,  15,  135,   90),
    lat_c = c(  45,   45,  45,   45,  45,   30,  45,   45,    -15,   0,  -15,   30 )
  )

  # Find the nearest window center
  centers_df$dist2 <- (centers_df$lon_c - lon)^2 + (centers_df$lat_c - lat)^2
  nearest <- centers_df[which.min(centers_df$dist2), ]

  # Define map extent (±45° from center)
  lon_min <- max(-180, nearest$lon_c - 45)
  lon_max <- min(180,  nearest$lon_c + 45)
  lat_min <- max(-90,  nearest$lat_c - 45)
  lat_max <- min(90,   nearest$lat_c + 45)

  # Load world map
  world <- ne_countries(scale = "medium", returnclass = "sf")

  # Create point geometry
  point_sf <- st_sf(
    data.frame(lon = lon, lat = lat),
    geometry = st_sfc(st_point(c(lon, lat)), crs = 4326)
  )

  # Plot
  p <- ggplot(data = world) +
    geom_sf(fill = "gray75", color = "gray60") +
    geom_sf(data = point_sf, color = point_color, size = point_size) +
    coord_sf(
      xlim = c(lon_min, lon_max),
      ylim = c(lat_min, lat_max),
      expand = FALSE
    ) +
    ggplot2::labs(title = title_text, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        size = text_size,
        vjust = 1
      ),
      axis.text.y = ggplot2::element_text(
        size = text_size,
        angle = 90,
        hjust = 0.5,
        vjust = 0.5,
        margin = margin(r = 2, l = 0)
      ),
      plot.title = ggplot2::element_text(size = text_size, face = "plain",
                                         margin = margin(b = 3)),
      plot.title.position = "panel",
      panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_line(color = "gray90", linewidth = 0.25)
    )

  return(p)
}
