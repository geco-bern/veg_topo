# ~ 6 min on UBELIX
# ----------------- Set up-----------------

library(sf)
library(terra)
library(furrr)
library(future)
library(ggplot2)
library(patchwork)
library(ggmap)

source(here::here("R/config.R"))
source(here::here("R/get_lonlat_extent.R")) # in create_aligned_template() to process lidar data
source(here::here("R/create_aligned_template.R"))
source(here::here("R/plot_dem.R"))
source(here::here("R/plot_vegh.R"))
source(here::here("R/plot_twi.R"))
source(here::here("R/plot_tri.R"))

source(here::here("R/plot_google_img.R"))
source(here::here("R/plot_cor_twi_vegh.R"))
source(here::here("R/plot_cor_pval.R"))
source(here::here("R/plot_hex_scatter.R"))
source(here::here("R/plot_r_H_R.R"))

# Set worker numbers for different system
hostname <- trimws(tolower(system("hostname", intern = TRUE)))
if (hostname == "dash") workers = 4 else workers = 12
message("→ using ", workers, " workers")

# All ALS_MAX data
file_bases <- c("30SWH","30TWN","30UVD","31UFT","32TMS",
                "32TMT", "34WFS", "35VMF","35VNL","35WMR","35WNT",
                "06VXR", "08WNA", "12VUN", "16PHS", "32MPC", "32MQE",
                "06WVT", "10TER", "13UEA", "16SGE", "32MPE", "32NNF")

# ----------------- Funtion to process single files-----------------
plot_validation_results <- function(file_base) {

  message("⭐️⭐️⭐️ Plotting: ", file_base, " ⭐️⭐️⭐️")
  # ============================================================================
  # 1. LOAD DATA AND CONFIGURATION
  # ============================================================================

  tryCatch({
    # ploting arguments
    text_size = 14
    x_step = 0.2
    y_step = 0.2
    res_out = 0.2

    # ----- Check if input file exists -----
    # LiDAR input path and alignment grid
    lidar_path <- file.path(lidar_asl_dir, paste0(file_base, ".tif"))
    lidar_path <- if (file.exists(lidar_path)) lidar_path else file.path(lidar_lvis_dir, paste0(file_base, ".tif"))
    if (!file.exists(lidar_path)) return(message("❌ LiDAR file not found"))

    # ----- Region info for ploting -----
    reg_extent <- ext(create_aligned_template(lidar_path, res_out = res_out))
    message("Extent:",reg_extent[1], ", ", reg_extent[2],", ",reg_extent[3],", ",reg_extent[4], " (xmin,xmax,ymin,ymax)")

    aspect_ratio <- (reg_extent[4] - reg_extent[3]) / (reg_extent[2] - reg_extent[1])
    message(sprintf("Aspect ratio (height/width): %.3f", aspect_ratio))

    # ----- Input file paths -----
    hlidar_file <- file.path(h_validation_dir, paste0(file_base, "_hlidar.tif"))
    hlang_file <- file.path(h_validation_dir, paste0(file_base, "_hlang.tif"))
    twi_file <- file.path(h_validation_dir, paste0(file_base, "_twi.tif"))
    dem_file <- file.path(h_validation_dir, paste0(file_base, "_dem.tif"))
    rin_file <- file.path(h_validation_dir, paste0(file_base, "_rin.tif"))

    # Check if required input files exist
    required_files <- c(hlidar_file, hlang_file, twi_file, dem_file, rin_file)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0) {
      message("❌ Missing required files:")
      for (f in missing_files) message("   - ", f)
      return(FALSE)
    }

    # ----- Output correlation files -----
    r_hlidar_twi_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlidar_twi.tif"))
    r_hlang_twi_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlang_twi.tif"))
    r_hlidar_rin_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlidar_rin.tif"))
    r_hlang_rin_file <- file.path(h_validation_dir, paste0(file_base, "_r_hlang_rin.tif"))

    # Check if correlation files exist
    cor_files <- c(r_hlidar_twi_file, r_hlang_twi_file, r_hlidar_rin_file, r_hlang_rin_file)
    missing_cor <- cor_files[!file.exists(cor_files)]
    if (length(missing_cor) > 0) {
      message("⚠️ Warning: Missing correlation files, output plots may be incomplete:")
      for (f in missing_cor) message("   - ", f)
    }

    # ============================================================================
    # 2. DEFINE BASE THEME
    # ============================================================================

    base_theme <- ggplot2::theme(
      aspect.ratio = aspect_ratio,

      legend.position = "right",
      panel.background = element_rect(fill = NA, color = NA),
      legend.background = element_rect(fill = NA, color = NA),
      legend.box.background = element_rect(fill = NA, color = NA),
      legend.text = ggplot2::element_text(size = text_size * 0.9,
                                          angle = 90,
                                          hjust = 0.5, vjust = 0.5,
                                          margin = margin(r = 0, b = 0, l = 0)),
      legend.title = ggplot2::element_text(size = text_size,
                                           angle = 0, hjust = 0, vjust = 1),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, 0, 0),

      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),

      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),

      axis.line = element_blank(),

      panel.spacing = unit(0, "cm"),
      panel.border = ggplot2::element_rect(linewidth = 0.3, fill = NA),
      plot.title = ggplot2::element_text(
        size = text_size * 1.2,
        face = "plain",
        margin = margin(b = 0)
      ),
      plot.title.position = "panel"
    )

    my_guides <- guides(fill = guide_colorbar(barwidth = 0.8, barheight = 10))

    # ============================================================================
    # 3. PLOT INPUT DATA
    # ============================================================================
    message("📊 Plotting input data...")

    p_google <- plot_google_img(extent = reg_extent, title_text = "   Google Satellite", text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme

    p_vegh_lang <- plot_vegh(hlang_file, extent = reg_extent,  title_text = expression("   30-m Lang" * italic(H)), text_size = text_size, x_step = x_step, y_step = y_step, color_limits = c(0, 70)) +
      base_theme + my_guides

    p_vegh_lidar <- plot_vegh(hlidar_file, extent = reg_extent, title_text = expression("   30-m Lidar" * italic(H)), text_size = text_size, x_step = x_step, y_step = y_step, color_limits = c(0, 70)) +
      base_theme + my_guides

    p_dem <- plot_dem(dem_file, extent = reg_extent, title_text = "   30-m Elevation", text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides

    p_twi <- plot_twi(twi_file, extent = reg_extent, title_text = "   30-m TWI", text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides

    p_rin <- plot_tri(rin_file, extent = reg_extent, title_text = "   30-m Rᵢₙ", text_size = text_size, x_step = x_step, y_step = y_step) +
      base_theme + my_guides

    # Combine input plots
    final_plot_input <- patchwork::wrap_plots(
      p_google, p_vegh_lang, p_vegh_lidar,
      p_dem, p_twi, p_rin,
      ncol = 3, nrow = 2
    ) +
      patchwork::plot_annotation(
        title = sprintf(
          "%s: Region: xmin=%s xmax=%s ymin=%s ymax=%s",
          file_base,
          reg_extent[1], reg_extent[2], reg_extent[3], reg_extent[4]
        ),
        tag_levels = "a"
      ) &
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank()
      )

    # Create figures directory if it doesn't exist
    fig_dir <- here::here("data/figures")
    if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

    # Save input plot
    out_file1 <- here::here("data/figures", paste0("6_01_validate_H_", file_base, "_input.png"))
    ggsave(filename = out_file1, plot = final_plot_input, width = 14, height = 8, dpi = 600)
    message("✅ Saved input plot: ", out_file1)

    # ============================================================================
    # 4. PLOT OUTPUT DATA (CORRELATIONS)
    # ============================================================================
    message("📈 Plotting correlation outputs...")

    # Only create correlation plots if files exist
    if (all(file.exists(c(r_hlidar_twi_file, r_hlang_twi_file,
                          r_hlidar_rin_file, r_hlang_rin_file)))) {

      p_rA_lidar <- plot_cor_twi_vegh(r_hlidar_twi_file, extent = reg_extent, title_text = bquote("   500-m " * r[.("Hlidar")*","*.("TWI")]), text_size = text_size,  x_step = x_step, y_step = y_step) +
        base_theme + my_guides

      p_rA_lang <- plot_cor_twi_vegh(r_hlang_twi_file, extent = reg_extent, title_text = bquote("   500-m " * r[.("Hlang")*","*.("TWI")]), text_size = text_size, x_step = x_step, y_step = y_step) +
        base_theme + my_guides

      p_rB_lidar <- plot_r_H_R(r_hlidar_rin_file, extent = reg_extent, title_text = bquote("   500-m " * r[.("Hlidar")*","*.("Rᵢₙ")]), text_size = text_size, x_step = x_step, y_step = y_step) +
        base_theme + my_guides

      p_rB_lang <- plot_r_H_R(r_hlang_rin_file, extent = reg_extent, title_text = bquote("   500-m " * r[.("Hlang")*","*.("Rᵢₙ")]), text_size = text_size, x_step = x_step, y_step = y_step) +
        base_theme + my_guides

      # Combine output plots
      final_plot_output <- patchwork::wrap_plots(
        p_google, p_rA_lang, p_rA_lidar,
        p_dem, p_rB_lang, p_rB_lidar,
        ncol = 3, nrow = 2
      ) +
        patchwork::plot_annotation(
          title = sprintf(
            "%s: Region: xmin=%s xmax=%s ymin=%s ymax=%s",
            file_base,
            reg_extent[1], reg_extent[2], reg_extent[3], reg_extent[4]
          ),
          tag_levels = "a"
        ) &
        theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_blank(),
          legend.background = element_blank(),
          legend.box.background = element_blank()
        )

      # Save output plot
      out_file2 <- here::here("data/figures", paste0("6_01_validate_H_", file_base, "_output.png"))
      ggsave(filename = out_file2, plot = final_plot_output, width = 14, height = 8, dpi = 600)
      message("✅ Saved output plot: ", out_file2)

    } else {
      message("⚠️ Skipping output plots: missing correlation files")
    }

    message("🎉 Plotting completed successfully for: ", file_base)
    return(TRUE)

  }, error = function(e) {
    message("\n❌ Error in plotting for ", file_base, ": ", e$message)
    return(FALSE)
  })
}


# ----------------- Parallel execution for all tiles-----------------
# Set up cluster plan
plan(cluster, workers = workers)

# Run in parallel
tictoc::tic("🚀 Parallel processing of tiles")

results <- future_map(
  file_bases,
  plot_validation_results,
  .progress = FALSE,
  .options = furrr_options(seed=TRUE)
)

plan(sequential)
tictoc::toc()

# Summarize results
success_count <- sum(unlist(results))
fail_count <- length(results) - success_count
message(sprintf("✅ Completed: %d succeeded, ❌ %d failed.", success_count, fail_count))

# ------ Single tile check (optional) ---------------------------------------
