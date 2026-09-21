# ~ 6.5 min on UBELIX
# ---------- SetUp -------------------------------------------------------------

library(terra)
library(ggplot2)
library(tidyterra)
library(scico)
library(rnaturalearth)
library(sf)

source(here::here("R/config.R"))
source(here::here("R/plot_rin.R"))

# ---------- Load Data ----------------------------------------------------------
# Load coastline vector data (for plotting reference only)
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# Load raster data
rin_450m <- terra::rast(rin_450m_mosaic_path)

rin_450m_agg <-  terra::aggregate(rin_450m, fact = 12)

# ------- Plot radiaiton index  ---------------------

p_te <- plot_rin(
  input = rin_450m_agg,
  extent = ext_global,
  text_size = 12,
  x_step = 30,
  y_step = 30
) +
  guides(fill = guide_colorbar(
  title.position = "left",
  barwidth = grid::unit(0.1, "in"),
  barheight = grid::unit(5, "in")
)) +
  geom_sf(data = coast,
          colour = 'black',
          linewidth = 0.1) +
  coord_sf(
    clip = "on"
  )

# save
ggsave(
  filename = file.path(here::here("data/figures/1_03_radiation_index_map.png")),
  plot = p_te,
  width = 14,
  height = 7,
  dpi = 600,
  units = "in"
)



