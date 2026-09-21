# ~84 min: UBELIX, 20 cores, 500G
# ~30 min: UBELIX, 35 cores, 1000G

# ------Load required libraries and custom functions----------------------------

library(terra)
library(dplyr)
library(furrr)

library(ggplot2)
library(tidyterra)
library(patchwork)
# library(ggmap)
library(khroma)
library(RColorBrewer)
library(rnaturalearth)
library(sf)

source(here::here("R/config.R"))
source(here::here("R/create_spatial_windows.R"))
source(here::here("R/calculate_correlation_bywin.R"))
source(here::here("R/mosaic_tiles.R"))
source(here::here("R/raster_preprocess_save.R"))
source(here::here("R/aggregate_topography.R"))
source(here::here("R/plot_dem.R"))
# source(here::here("R/plot_aspect.R"))
# source(here::here("R/plot_slope.R"))
source(here::here("R/plot_vegh.R"))
source(here::here("R/plot_twi.R"))
source(here::here("R/plot_r_H_R.R"))
source(here::here("R/plot_cor_twi_vegh.R"))
# source(here::here("R/plot_cor_pval.R"))
# source(here::here("R/plot_sw_in.R"))
source(here::here("R/plot_rin.R"))
# source(here::here("R/plot_hex_scatter.R"))
source(here::here("R/plot_single_sample_location.R"))
# source(here::here("R/plot_google_img.R"))
source(here::here("R/plot_fused.R"))
source(here::here("R/plot_kg_class.R"))


# Set workers, but need `source(here::here("R/config.R"))` first
if (hostname == "dash") workers = 2 else workers = 35
message("→ using ", workers, " workers")
# ------ File Configuration ---------------------------------------------

if (!dir.exists(r_H_R_tiles_dir)) {
  dir.create(r_H_R_tiles_dir, recursive = TRUE)
  message("Directory created: ", r_H_R_tiles_dir)
}

# ------------ function to process single tile----------------------------------
process_r_H_R_5km <- function(tile_row, output_dir = r_H_R_tiles_dir,
                              text_size = 12) {

  tryCatch({
    # --- Tile info ---
    tile_id <- tile_row$tile_id
    tile_extent <- terra::ext(tile_row$xmin, tile_row$xmax, tile_row$ymin, tile_row$ymax)
    tile_xmid <- (tile_row$xmin + tile_row$xmax)/2
    tile_ymid <- (tile_row$ymin + tile_row$ymax)/2

    x_step <- tile_row$xmax - tile_row$xmin
    y_step <-  tile_row$ymax - tile_row$ymin

    tictoc::tic(paste0("Processing tile: ", tile_id))
    t0 <- Sys.time()

    # --- TWI Raster ---
    twi_rc <- terra::rast(twi_450m_mosaic_clean_path) |> terra::crop(tile_extent)
    twi_tmp_path <- file.path(tempdir(), paste0("tile_", tile_id, "_twi_450m.nc"))
    terra::writeCDF(twi_rc, twi_tmp_path, varnames = "twi", overwrite = TRUE)
    rm(twi_rc); gc()
    twi_rc <- terra::rast(twi_tmp_path)
    names(twi_rc) <- "twi"
    message("Saved temporary TWI raster: ", twi_tmp_path)

    # --- Vegetation Height Raster ---
    vegh_rc <- terra::rast(vegh_450m_mosaic_path) |> terra::crop(tile_extent)
    vegh_tmp_path <- file.path(tempdir(), paste0("tile_", tile_id, "_vegh_450m.nc"))
    terra::writeCDF(vegh_rc, vegh_tmp_path, varnames = "vegh",overwrite = TRUE)
    rm(vegh_rc); gc()
    vegh_rc <- terra::rast(vegh_tmp_path)
    names(vegh_rc) <- "vegh"
    message("Saved temporary vegetation height raster: ", vegh_tmp_path)

    # --- Radiation index of terrain effect ---
    rin_rc <- terra::rast(rin_450m_mosaic_path) |> terra::crop(tile_extent)
    terrain_tmp_path <- file.path(tempdir(), paste0("tile_", tile_id, "_rin_450m.nc"))
    terra::writeCDF(rin_rc, terrain_tmp_path, varnames = "rin", overwrite = TRUE)
    rm(rin_rc); gc()
    rin_rc <- terra::rast(terrain_tmp_path)
    names(rin_rc) <- "rin"
    message("Saved temporary radiation index raster: ", terrain_tmp_path)

    # --- Stack and correlation ---
    stacked <- c(rin_rc, vegh_rc)
    df_win <- create_spatial_windows(stacked, value_vars = c("rin", "vegh"), dwin = 0.05)
    df_cor <- calculate_correlation_bywin(df_win, x = "rin", y = "vegh")

    # --- Save correlation ---
    cor_r <- terra::rast(df_cor[, c("lon_mid", "lat_mid", "correlation")], type="xyz", crs="EPSG:4326")
    names(cor_r) <- "r_H_R"
    cor_nc_path <- file.path(output_dir, paste0("tile_", tile_id, "_r_H_R_5km_map.nc"))
    terra::writeCDF(cor_r, cor_nc_path, overwrite = TRUE)
    message("Saved: ", cor_nc_path)

    # --- Save p-value ---
    pval_r <- terra::rast(df_cor[, c("lon_mid", "lat_mid", "cor_pval")], type="xyz", crs="EPSG:4326")
    names(pval_r) <- "pval_r_H_R"
    pval_nc_path <- file.path(output_dir, paste0("tile_", tile_id, "_r_H_R_5km_pval.nc"))
    terra::writeCDF(pval_r, pval_nc_path, overwrite = TRUE)
    message("Saved: ", pval_nc_path)

    # --- reset theme for plots ---
    re_theme <- list(
      guides(fill = guide_colorbar(barwidth = 0.8, barheight = 12)),
      ggplot2::coord_sf(
        clip = "on"
      ),
      ggplot2::theme(
        aspect.ratio = 1,
        legend.position = "right",
        legend.text = ggplot2::element_text(size = text_size * 0.9, angle = 90,
                                            hjust = 0.5, vjust = 0.5,
                                            margin = margin(r = 0, l = 2)),
        legend.title = ggplot2::element_text(size = text_size,
                                             angle = 90,   hjust = 0, vjust = 0.5 ),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 2, 0, -8),
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
          margin = margin(r = 2, l = -15) # important set to reduce the space between 2 plot
        ),
        panel.spacing = unit(0, "cm"),
        panel.border = ggplot2::element_rect(linewidth = 0.3, fill = NA),
        plot.margin = margin(0, 0, 0, 0),
        plot.title = ggplot2::element_text(
          size = text_size * 1.2,
          face = "plain",
          margin = margin(b = 0)
        ),
        plot.title.position = "panel"
      )
    )

    # The left indentation of the y-axis label should be removed from the first image in each row;
    # otherwise, the y-axis text will be clipped.
    re_theme_left <- ggplot2::theme(axis.text.y = ggplot2::element_text(
      margin = margin(r = 2, l = 0) # important set to left plot to show full y_text
    ))

    # ---- Generate plots ----
    p_dem <- plot_dem(dem_450m_mosaic_path, extent = tile_extent, title_text = "450 m: Elevation", text_size = text_size, x_step = x_step, y_step = y_step)  + re_theme
    p_vegh <- plot_vegh(vegh_rc, extent = tile_extent, title_text = expression("450 m: " * italic(H)[veg]), text_size = text_size, x_step = x_step, y_step = y_step)  + re_theme
    p_twi <- plot_twi(twi_rc, extent = tile_extent, title_text = "450 m: TWI", text_size = text_size, x_step = x_step, y_step = y_step)  + re_theme + re_theme_left

    p_rin <- plot_rin(rin_rc, extent = tile_extent, title_text = "450 m: Radiation index",  text_size = text_size, x_step = x_step, y_step = y_step) + re_theme
    p_rA <- plot_cor_twi_vegh(cor_twi_vegh_mosaic_file, extent = tile_extent,  title_text <- bquote("5 km: Pearson's " * r[.("H")*","*.("TWI")]), text_size = text_size, x_step = x_step, y_step = y_step)  + re_theme + re_theme_left
    p_rB <- plot_r_H_R(cor_nc_path, extent = tile_extent, title_text = bquote("5 km: Pearson's " * r[.("H")*","*.("Rᵢₙ")]), text_size = text_size, x_step = x_step, y_step = y_step) + re_theme

    p_fused <- plot_fused(fused_5km_file, extent = tile_extent, text_size = text_size, x_step = x_step, y_step = y_step) + re_theme
    p_kg <- plot_kg_class(kg_present_0p083_file, kg_legend_file, extent = tile_extent, text_size = text_size, x_step = x_step, y_step = y_step) + ggplot2::theme(aspect.ratio = 1)
    p_location <- plot_single_sample_location(tile_xmid, tile_ymid,  tile_id, text_size = text_size)  + ggplot2::theme(aspect.ratio = 1)

    # p_slope <- plot_slope(slope_450m_mosaic_path, extent = tile_extent, text_size = text_size) + re_theme
    # p_aspect <- plot_aspect(aspect_450m_mosaic_path, extent = tile_extent, text_size = text_size) + re_theme
    # p_p <- plot_cor_pval(pval_r, extent = tile_extent,  title_text = "5km: Pearson's p value (H～Rᵢₙ)", text_size = text_size) + re_theme

    # p_google <- plot_google_img(extent = tile_extent) + ggplot2::theme(aspect.ratio = 1)
    # p_scatter <- plot_hex_scatter(df_win, x_var = "rin", y_var = "vegh", x_text = "Radiation index", y_text = "Vegetation height (m)", text_size = text_size) + ggplot2::theme(aspect.ratio = 1)

    # ---- Combine plots ----
    final_plot <- (
      (p_twi + p_dem + p_rin) /
        (p_rA + p_vegh + p_rB)  /
        (p_kg + p_fused + p_location)) +
      plot_annotation(title = tile_id) +
      plot_layout(heights = c( 1, 1, 1))

    # ---- Save plot ----
    out_file <- here::here(file.path(paste0("data/figures/3_01_tile_", tile_id, "_H_Rin_plots.png")))
    ggsave(filename = out_file, plot = final_plot, width = 14, height = 14.2, dpi = 600)

    # ---- Memory cleanup ----
    # list all objects need to be remove
    rm(
      twi_rc, vegh_rc, rin_rc,
      stacked, df_win, df_cor, cor_r,pval_r,
      p_dem, p_vegh, p_rA, p_rB,
      p_kg, p_fused, p_location,
      final_plot)
    gc(verbose = FALSE)

    # --- Print proccessed time ---
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("Region %s completed [%.1f mins]", tile_id, elapsed_mins))
    tictoc::toc()

    return(TRUE)
  }, error = function(e) {
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("❌ Tile %s failed after %.1f mins: %s", tile_row$tile_id, elapsed_mins, e$message))
    return(FALSE)
  })
}


# ----------------- Parallel execution -----------------
# Set up cluster plan
plan(cluster, workers = workers)

# load tiles info
tiles_info <- readRDS(valid_tiles_info_path)

# Run in parallel
tictoc::tic("🚀 Parallel processing of tiles")
results <- future_map(seq_len(nrow(tiles_info)),
                      function(i) process_r_H_R_5km(tiles_info[i, ]),
                      .progress = FALSE,
                      .options = furrr_options(seed=TRUE))
plan(sequential)
tictoc::toc()

# Summarize results
success_count <- sum(unlist(results))
fail_count <- length(results) - success_count
message(sprintf("✅ Completed: %d succeeded, ❌ %d failed.", success_count, fail_count))

# -------- Combination ---------------------------------------------------------
# mosacing the r(H~Rin) map
mosaic_tiles(
  input_dir   = r_H_R_tiles_dir,
  output_file = r_H_R_5km_path,
  pattern = "*_map.nc",
  varname = "r_H_R")

# mosacing the pval ofr(H~Rin) map
mosaic_tiles(
  input_dir   = r_H_R_tiles_dir,
  output_file = pval_r_H_R_5km_path,
  pattern = "*_pval.nc",
  varname = "pval_r_H_R")

# ---------- Delete intermediate data ------------------------------------------
# List all files in the directory r_H_R_tiles_dir that match "*.nc"
# If there are any files found, delete them

# cor_5km_tiles_path <- fs::dir_ls(path = r_H_R_tiles_dir, glob = "*.nc")
# if (length(cor_5km_tiles_path) > 0) file.remove(cor_5km_tiles_path)

# # ------ Single tile test ---------------------------------------------
# # test for 1 tile
# tiles_info <- readRDS(valid_tiles_info_path)
# process_r_H_R_5km(tiles_info[40,])

# # ------ Smaller region  test ---------------------------------------------
# regA_info <- data.frame(
#   tile_id = c("Aletsch_glacier"),
#   ymin = c(46.8),
#   ymax = c(47),
#   xmin = c(7.8),
#   xmax = c(8),
#   sample_id = c(1)
# )
#
# # center location
# regA_info$xmid <- (regA_info$xmax + regA_info$xmin) / 2
# regA_info$ymid <- (regA_info$ymax + regA_info$ymin) / 2
#
# tile_row <- regA_info[1, ]
#
# output_dir = r_H_R_tiles_dir
# text_size = 12
# x_step = 0.05
# y_step = 0.05



