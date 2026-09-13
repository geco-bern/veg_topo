
# ~ UBELIX: 8.5 min for 34 sample regions
# ------ Load required libraries -------------------------------------------------------------
library(terra)
library(dplyr)
library(tidyr)
library(fs)
library(stringr)
library(parallel)
library(future)
library(furrr)
library(meteoland)

library(ggplot2)
library(patchwork)
library(ggmap)
library(khroma)
library(RColorBrewer)
library(rnaturalearth)
library(sf)

# ------ Load helper and custom functions -----------------------------------------------------
source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R")) # with source(here::here("R/save_raster.R"))
source(here::here("R/df_to_raster.R")) # with source(here::here("R/save_raster.R"))
source(here::here("R/create_aligned_template.R"))
source(here::here("R/calc_meteoland_sw_in.R"))
source(here::here("R/calc_tri.R"))
source(here::here("R/mosaic_tiles.R"))
source(here::here("R/extent_to_tile_ids.R"))
source(here::here("R/create_spatial_windows.R"))
source(here::here("R/calculate_correlation_bywin.R"))

source(here::here("R/plot_dem.R"))
source(here::here("R/plot_vegh.R"))
source(here::here("R/plot_twi.R"))
source(here::here("R/plot_tri.R"))

source(here::here("R/plot_google_img.R"))
source(here::here("R/plot_cor_twi_vegh.R"))
source(here::here("R/plot_cor_pval.R"))
source(here::here("R/plot_hex_scatter.R"))
source(here::here("R/plot_r_H_R.R"))

source(here::here("R/plot_single_sample_location.R"))
source(here::here("R/plot_kg_class.R"))
source(here::here("R/plot_fused.R"))
source(here::here("R/plot_kmeans_map.R"))
source(here::here("R/plot_mi.R"))
source(here::here("R/plot_cci_land_cover.R"))
source(here::here("R/plot_scatter_r_validation.R"))
# Optional (uncomment if needed)

source(here::here("R/plot_aspect.R"))
source(here::here("R/plot_slope.R"))
source(here::here("R/plot_sw_in.R"))

# ---------------  Set worker numbers for different system ---------------------
# Set worker configuration based on system
if (hostname == "dash") {
  chunk_size <- 1000
  workers <- 16
} else {
  chunk_size <- 10000
  workers <- 100
}
message("→ Using ", workers, " workers and chunk_size = ", chunk_size, " for parallel processing.\n")

# --------------- File Configuration -------------------------------------------
if (!dir.exists(reg_validate_dir)) {
  dir.create(reg_validate_dir, recursive = TRUE)
  message("Directory created: ", reg_validate_dir)
}

# --------------- Main Processing Function -------------------------------------
process_reg_500m <- function(reg_row, output_dir = reg_validate_dir, text_size = 14) {
  t0 <- Sys.time()

  # ============ Setup per regions ============================================

  # --- Region info ---
  reg_id <- reg_row$strata_label
  reg_extent <- terra::ext(reg_row$xmin, reg_row$xmax, reg_row$ymin, reg_row$ymax)
  reg_xmid <- (reg_row$xmin +reg_row$xmax) / 2
  reg_ymid <- (reg_row$ymin +reg_row$ymax) / 2
  grid_30m <- create_aligned_template(reg_extent, res_out = 0.00025)

  # ----- Plot config -----
  x_step <- reg_extent[2] - reg_extent[1]
  y_step <-  reg_extent[4] - reg_extent[3]
  aspect_ratio <- (reg_extent[4] - reg_extent[3]) / (reg_extent[2] - reg_extent[1])

  # ----- Output file paths -----
  vegh_file <- file.path(output_dir, paste0(reg_id, "_vegh_30m.tif"))
  twi_file <- file.path(output_dir, paste0(reg_id, "_twi_30m.tif"))
  dem_file <- file.path(output_dir, paste0(reg_id, "_dem_30m.tif"))
  slope_file <- file.path(output_dir, paste0(reg_id, "_slope_degrees_30m.tif"))
  aspect_file <- file.path(output_dir, paste0(reg_id, "_aspect_degrees_30m.tif"))
  tri_file <- file.path(output_dir, paste0(reg_id, "_tri_30m.tif"))

  r_h_twi_file <- file.path(output_dir, paste0(reg_id, "_r_h_twi_30m_500m.tif"))
  p_h_twi_file <- file.path(output_dir, paste0(reg_id, "_p_h_twi_30m_500m.tif"))
  r_h_tri_file <- file.path(output_dir, paste0(reg_id, "_r_h_tri_30m_500m.tif"))
  p_h_tri_file <- file.path(output_dir, paste0(reg_id, "_p_h_tri_30m_500m.tif"))

  r_h_twi_30_5000_file <- file.path(output_dir, paste0(reg_id, "_r_h_twi_30_5000.tif"))
  p_h_twi_30_5000_file <- file.path(output_dir, paste0(reg_id, "_p_h_twi_30_5000.tif"))
  r_h_tri_30_5000_file <- file.path(output_dir, paste0(reg_id, "_r_h_tri_30_5000.tif"))
  p_h_tri_30_5000_file <- file.path(output_dir, paste0(reg_id, "_p_h_tri_30_5000.tif"))

  tryCatch({

    # ============ Start Processing ============================================
    tictoc::tic(paste0("⭐️⭐️⭐️ Processing: ", reg_id, " ⭐️⭐️⭐"))
    t0 <- Sys.time()

    # --- TWI Raster ---
    twi_rc <- terra::rast(twi_30m_path) |> terra::crop(reg_extent)
    twi_rc <- twi_rc / 100
    save_raster(twi_rc, twi_file)

    # --- Vegetation Height Raster ---
    vegh_rc <- extent_to_tile_ids(reg_extent, tile_size = 3, return_raster = TRUE,
                                  source = "lang_vegh_10m", tiles_dir = vegh_10m_tiles_dir)

    # Set 0 as NA value (0m canopy height represents not vegetated or water according to Lang et al. (2019))
    # Aggregate and resample using TWI data from Ho et al. (2025)
    vegh_rc <- raster_preprocess_save(
      input = vegh_rc,
      output = vegh_file,
      target = grid_30m,
      na_value = 0,
      fun = mean,
      varname = "vegh",
      if_aggregate = TRUE,
      if_resample = TRUE,
      if_return_raster = TRUE
    )

    # --- Stack and correlation using 30-m input data and 500m window---
    stacked <- c(twi_rc, vegh_rc)
    df_win <- create_spatial_windows(stacked, value_vars = c("twi", "vegh"), dwin = 0.005)
    df_cor <- calculate_correlation_bywin(df_win, x = "twi", y = "vegh")

    # --- Save correlation ---
    r_h_twi <- df_to_raster(df_cor, "lon_mid", "lat_mid", "correlation", grid_30m) |>
      save_raster(r_h_twi_file)
    p_h_twi <- df_to_raster(df_cor, "lon_mid", "lat_mid", "cor_pval", grid_30m) |>
      save_raster(p_h_twi_file)

    # --- Stack and correlation using 30-m input data and 5km window---
    df_win <- create_spatial_windows(stacked, value_vars = c("twi", "vegh"), dwin = 0.05)
    df_cor <- calculate_correlation_bywin(df_win, x = "twi", y = "vegh")

    r_h_twi_30_5000 <- df_to_raster(df_cor, "lon_mid", "lat_mid", "correlation", grid_30m) |>
      save_raster(r_h_twi_30_5000_file)

    p_h_twi_30_5000 <- df_to_raster(df_cor, "lon_mid", "lat_mid", "cor_pval", grid_30m) |>
      save_raster(p_h_twi_30_5000_file)

    # --- Elevation, slope and aspect ---
    # get raw dem and then calculate slope, aspect and TRI
    dem_ex <- extent_to_tile_ids(reg_extent, tile_size = 1, return_raster = TRUE,
                                 source = "COP30", tiles_dir = COP30_dir, if_crop = FALSE)
    topo_vars <- calc_tri(dem_ex, grid_30m, chunk_size, workers)

    # Save DEM, Slope, Aspect
    save_raster(topo_vars$dem, dem_file)
    save_raster(topo_vars$slope, slope_file)
    save_raster(topo_vars$aspect, aspect_file)
    save_raster(topo_vars$tri, tri_file)

    # --- Stack and Correalation using 30-m input data and 0.005-degree window ---
    stacked <- c(topo_vars$tri, vegh_rc)
    names(stacked) <- c("tri", "vegh")
    df_win <- create_spatial_windows(stacked, value_vars = c("tri", "vegh"), dwin = 0.005)
    df_cor <- calculate_correlation_bywin(df_win, x = "tri", y = "vegh")

    # --- Save correlation ---
    r_h_tri <- df_to_raster(df_cor, "lon_mid", "lat_mid", "correlation", grid_30m) |>
      save_raster(r_h_tri_file)
    p_h_tri <- df_to_raster(df_cor, "lon_mid", "lat_mid", "cor_pval", grid_30m) |>
      save_raster(p_h_tri_file)

    # --- Stack and Correalation using 30-m input data and 0.05-degree window ---
    df_win <- create_spatial_windows(stacked, value_vars = c("tri", "vegh"), dwin = 0.05)
    df_cor <- calculate_correlation_bywin(df_win, x = "tri", y = "vegh")

    # --- Save correlation ---
    r_h_tri_30_5000 <- df_to_raster(df_cor, "lon_mid", "lat_mid", "correlation", grid_30m) |>
      save_raster(r_h_tri_30_5000_file)
    p_h_tri_30_5000 <- df_to_raster(df_cor, "lon_mid", "lat_mid", "cor_pval", grid_30m) |>
      save_raster(p_h_tri_30_5000_file)

    # --- reset theme for plots, change the layout for better visualization ---
    base_theme <- ggplot2::theme(
      aspect.ratio = aspect_ratio,

      legend.position = "right",
      legend.justification = c(0, 0), # 左对齐，垂直底部对齐
      panel.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      legend.box.background = element_rect(fill = NA, color = NA),
      legend.text = ggplot2::element_text(size = text_size * 0.9,
                                          angle = 0,
                                          hjust = 0, vjust = 0.5,
                                          margin = margin(b = 0) ),
      legend.title = ggplot2::element_text(size = text_size,
                                           angle = 0, hjust = 0, vjust = 1),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, -8),

      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      # axis.ticks.x = element_line(),
      axis.ticks.x = element_blank(),

      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      # axis.ticks.y =  element_line(),
      axis.ticks.y =  element_blank(),

      axis.line = element_blank(),

      panel.spacing = unit(0, "cm"),
      # panel.border = ggplot2::element_rect(linewidth = 0.3, fill = NA),
      plot.title = ggplot2::element_text(
        size = text_size * 1.2,
        face = "plain",
        margin = margin(b = 0)
      ),
      plot.title.position = "panel"
    )

    my_guides <- guides(fill = guide_colorbar(label.position = "right",  # 标签在右侧
                                              label.hjust = 0,
                                              barwidth = 0.8, barheight = 12,
                                              label.theme = element_text(
                                                size = text_size * 0.9,
                                                margin = margin(r = 2, l = 2)  # 统一右侧固定边距
                                              )))

    # re_theme_margin_r for the first fig in each rows
    # re_theme_margin_r <- ggplot2::theme( plot.margin = margin(t = 0, r = 20, b = 0, l = 0))

    # --- Plotting ---
    p_location <- plot_single_sample_location(reg_xmid, reg_ymid, title_text = "    Location", text_size = text_size) +
      ggplot2::theme(aspect.ratio = aspect_ratio, plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_google <- plot_google_img(extent = reg_extent, title_text = "    Google Satellite", text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + ggplot2::theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_dem <- plot_dem(dem_file, extent = reg_extent, title_text = "    30-m Elevation", text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides + ggplot2::theme( plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_vegh <- plot_vegh(vegh_file, extent = reg_extent, title_text = expression("    30-m " * italic(H)), text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides + ggplot2::theme( plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_twi <- plot_twi(twi_file, extent = reg_extent, title_text = "    30-m TWI", text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides + ggplot2::theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_rA <- plot_cor_twi_vegh(r_h_twi, extent = reg_extent,  title_text <- bquote("    500-m Pearson's " * r[.("H")*","*.("TWI")]), text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides + ggplot2::theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_tri <- plot_tri(tri_file, extent = reg_extent, title_text = "    30-m TRI",  text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides + ggplot2::theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_rB <- plot_r_H_R(r_h_tri, extent = reg_extent, title_text = bquote("    500-m Pearson's " * r[.("H")*","*.("TRI")]), text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides + ggplot2::theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    # ---- x=r30; y=r450 ----
    p_validA <- plot_scatter_r_validation(input_x = r_h_twi_30_5000_file,
                                          input_y = cor_twi_vegh_mosaic_file,
                                          title_text = bquote("    Comparison " * r[.("H")*","*.("TWI")]),
                                          text_size = text_size) +
      ggplot2::theme(legend.position = "none", plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_rA2 <- plot_cor_twi_vegh(cor_twi_vegh_mosaic_file, extent = reg_extent, title_text = bquote("    5-km Pearson's " * r[.("H")*","*.("TWI")]), text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides + ggplot2::theme( plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_validB <- plot_scatter_r_validation(r_h_tri_30_5000_file, r_H_R_5km_path, title_text = bquote("    Comparison " * r[.("H")*","*.("TRI")]), text_size = text_size) +
      ggplot2::theme(legend.position = "none", plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    p_rB2 <- plot_r_H_R(r_H_R_5km_path, extent = reg_extent, title_text = bquote("    5-km Pearson's " * r[.("H")*","*.("TRI")]), text_size = text_size, x_step = x_step, y_step = y_step)+
      base_theme + my_guides + ggplot2::theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

    final_plot1 <- patchwork::wrap_plots(
      p_location, p_google,  p_dem, p_vegh,
      p_validA, p_rA, p_rA2, p_twi,
      p_validB, p_rB, p_rB2, p_tri,
      ncol = 4, nrow = 3
    ) +
      patchwork::plot_annotation(
        title = sprintf(
          "%s:\n xmin = %s, xmax = %s, ymin = %s, ymax = %s",
          reg_id,
          reg_extent[1], reg_extent[2], reg_extent[3], reg_extent[4]
        ),
        tag_levels = "a",
        tag_suffix = ")"
      ) &
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        plot.tag.position = c(0.02, 0.97),
        plot.tag = element_text(size = text_size, face = "bold"),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
      )

    out_file1 <- here::here(file.path(paste0("data/figures/5_02_1_validate_", reg_id, "_12plots.png")))
    ggsave(filename = out_file1, plot = final_plot1, width = 14, height = 12.1, dpi = 600)


    # --- Print proccessed time ---
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("Region %s completed [%.1f mins]", reg_id, elapsed_mins))
    tictoc::toc()
    return(TRUE)

  }, error = function(e) {
    reg_id <- paste0(reg_row$strata_label, "_", reg_row$strata_label)
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("❌ Tile %s failed after %.1f mins: %s", reg_row$strata_label, elapsed_mins, e$message))
    return(FALSE)
  })
}

# ------------ validation region define ----------------------------------------
# all samples regions from 5_01
reg_info_all_samples <- readRDS(reg_sample_info_path) |>
  select(ends_with("label"), ends_with("min"), ends_with("max"), -starts_with("dem"))
print(reg_info_all_samples )

# # select 6 samples from 5_01 combined all MI bins and abs_lat bins

reg_info <- tribble(
  ~strata_label,                        ~ymin,   ~ymax,   ~xmin,     ~xmax,
  "Loetschental",                      46.4,    46.5,      7.8,       7.9,

  "b29_Amazonia_NE",                     -0.5,     0.5,     -55.5,       -54.5,
  "b27_Argentine",                        -36.5,  -35.5,     -60.5,      -59.5,
  "b30_Andes",                            -1,         0,     -78,      -77,

  "b28_Amazonia_NW",                         0,     1,     -62.5,       -61.5,
  "3d_waterlogged_pantanal",            -17.5,    -16.5,   -57.5,     -56.5,

  "hyper_arid_low_lat_rugged_relief",      -9,  -8.5,      38,      38.5,
  "b3_equatorial_rainforest_CongoBasin",  -1.0,    -0.5,     17.0,      17.5,
  "a1_waterlogged_pantanal",              -17.5,   -16.5,-57.5,     -56.5,
  "b25_Finland",                          67.5,    68.5,     25,        26,
  "b26_Monte_Alen_Guiana",                 1.2,       1.7,      9.8,        10.3,
  "a19_arctic_tundra_alaska",              68.5,    69,   -146,    -145.5,

  # ------------------ Subsurface flow validation (Fan et al., 2019) -----------
  "3a_desert_riparian",                 31,       32,      -110.5,   -109.5,
  "3b_mediterranean_california",        37,       38,      -122.5,   -121.5,
  "3c_forest_savanah_east_congo",       -5.5,     -4.5,     28,        29,

  "3e_white_sand_amazon",               3,         4,      -67.7,     -66.7,
  "3f_cool_wet_denmark",                56.7,     57.0,     10.0,      10.5,

  # ------------------ Additional validation from global results -----------

  "a2_amazon_floodplain_colombia",       0.5,     1.5,    -70,       -69,
  "a3_peru_western_amazon_terra_firme",  0,       1,      -71.5,     -70.5,
  "a4_forest_savanna_congo_mayombe",    -5.2,    -4.2,     17,        18,
  "a5_miombo_savanna_angola",           -8,      -7,       16.3,      17.3,
  "a6_zambezi_riparian_zambia",        -12,     -11,       18.5,      19.5,
  "a7_madagascar_highland_mosaic",     -21,     -20,       46.5,      47.5,
  "a8_hindu_kush_steppe_afghanistan",   34,      35,       66.5,      67.5,
  "a9_korean_temperate_forest",         35.7,    36.7,    127.6,     128.6,
  "a10_siberian_taiga_russia",          61,      62,      142.7,     143.7,
  "a11_borneo_lowland_rainforest",       0.2,     1.2,    112.8,     113.8,
  "a12_great_plains_grassland_usa",     38.5,    39.5,    -99,       -98,
  "a13_swiss_plateau_mixed_forest",     46.8,    47.8,      7,         8,
  "a14_french_alps_montane_forest",     44.5,    45.5,      6,         7,
  "a15_southwest_france_temperate",     43.6,    44.6,      2,         3,
  "a16_basque_coastal_forest_spain",    41.5,    42.5,     -2,        -1,
  "a17_yungas_cloud_forest_bolivia",   -23.3,   -22.3,    -65.3,     -64.3,
  "a18_gran_chaco_dry_savanna",        -26.6,   -25.6,    -60,       -59,
  "a20_jiangxi_lowlands_china",         28,      29,      114,       115,
  "a21_fujian_jiangxi_hills_china",     25.5,    26.5,    115.3,     116.3,

  # ------------------ Slope aspect validation (Fan et al., 2019) --------------
  # --- North America (Arid & Temperate regions) ---
  "2a_arid_western_Texas",              31.5,    32.5,   -105.0,    -104.5,
  "2b_seasonal_arid_California_foothills",37.5,  38.5,   -120.5,    -119.5,
  "2c_seasonal_arid_Idaho_Basin",       43.0,    44.0,   -115.8,    -114.8,
  "2d_winter_snowpacks_Mount_Cramer",   43.8,    44.8,   -115.5,    -114.5,

  # --- Northern high-latitude regions ---
  "2e_high_latitude_southwest_Yukon",   61.5,    62.5,   -138.5,    -137.5,
  "2f_boreal_forest_Alaska",            64.0,    65.0,   -148.0,    -147.0,

  # ------------------ Additional validation from global results ----------------

  "b2_Parc_National_Suisse",            46.5,    46.8,     10.1,      10.4,

  "b4_WesternAlps_MontBlanc_massif",    45.0,    45.5,      6.5,       7.0,
  "b5_subtropical_Yunnan_region",       35.5,    36.0,     97.0,      97.5,
  "b6_Orinoco_transition_forest",       -2.5,    -1.5,    120,       121,
  "b7_borneo_lowland_rainforest",        1,       2,      112.8,     113.8,
  "b8_victoria_coastal_forest_australia",-38,    -37,     146.2,     147.2,
  "b9_equatorial_forest_congo",         -0.5,     0.5,     16,        17,
  "b10_lesotho_mountain_grassland",    -30,     -29,      27.5,      28.5,
  "b11_ethiopian_highlands",             7.5,     8.5,     39,        40,
  "b12_bolivia_andes_cloud_forest",    -18.8,   -17.8,    -66.1,     -65.1,
  "b13_pantanal_wetland_brazil",       -19,     -18,     -57,       -56,
  "b14_patagonian_steppe_chile",       -47,     -46,     -69,       -68,
  "b15_altiplano_bolivia_peru",        -18,     -17,     -69.5,     -68.5,
  "b16_amazon_rainforest_brazil",       -1.5,    -0.5,    -62.7,     -61.7,
  "b17_arctic_tundra_alaska",           68.5,    69.5,   -146,      -145,
  "b18_boreal_forest_yukon",            62.6,    63.6,   -142.4,    -141.4,
  "b19_boreal_forest_nwt_canada",       63.5,    64.5,   -130.5,    -129.5,
  "b20_mississippi_river_valley_usa",   38.5,    39.5,    -89,       -88,
  "b21_western_anatolia_turkey",        37,      38,      28,        29,
  "b22_khorat_plateau_thailand",        15.5,    16.5,    103,       104,
  "b23_nepal_middle_hills",             27.2,    28.2,     81.8,      82.8,
  "b24_tianshan_mountains_xinjiang",    42.5,    43.5,     86,        87
)


# ----------- Process sample regionsand other validation regions-----------------------------
process_reg_500m(reg_info[1, ]) # Loetschental, smallest one, best for single region test

# for (i in 1:6) {
#   process_reg_500m(reg_info[i, ])ß
# }
#
# for (i in 7:13) {
#   process_reg_500m(reg_info[i, ])
# }
#
# for (i in 14:66) {
#   process_reg_500m(reg_info[i, ])
# }
#
# for (i in seq_len(nrow(reg_info_all_samples))) {
#   process_reg_500m(reg_info_all_samples[i, ]) # ~ 8h 43 min
# }
#
# for (i in seq_len(nrow(reg_info))) {
#   process_reg_500m(reg_info[i, ]) # ~ 8h 43 min
# }
