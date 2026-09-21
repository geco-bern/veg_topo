# ---- x=r30,y=r450 -----
plot_scatter_r_validation <- function(input_x,
                                      input_y,
                                      title_text = "Comparison r",
                                      x_text = "r(30)",
                                      y_text = "r(450)",
                                      text_size = 6
) {

  # ---- Load raster ----
  if (is.character(input_x)) input_x <- terra::rast(input_x)
  if (!inherits(input_x, "SpatRaster")) stop("input_x must be a SpatRaster or valid file path.")

  if (is.character(input_y)) input_y <- terra::rast(input_y)
  if (!inherits(input_y, "SpatRaster")) stop("input_y must be a SpatRaster or valid file path.")

  # ---- Handle extent and resampling ----
  input_y <- terra::resample(input_y, input_x, method = "near")

  # ---- Create dataframe ----
  stacked <- c(input_x, input_y)
  df <- as.data.frame(stacked, xy = FALSE, na.rm = FALSE)
  colnames(df) <- c("input_x", "input_y")
  df_clean <- na.omit(df)  # 移除NA值

  df_clean
  # ---- Create plot ----
  p <- ggplot2::ggplot(df_clean) +
    # 1:1 reference line (red dashed)
    ggplot2::geom_abline(intercept = 0, slope = 1,
                         color = "firebrick",
                         linetype = "dashed",
                         linewidth = 0.5) +
    # x=0 和 y=0 线
    ggplot2::geom_vline(xintercept = 0,
                        color = "gray80",
                        linetype = "solid",
                        linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0,
                        color = "gray80",
                        linetype = "solid",
                        linewidth = 0.8) +

    # Scatter points
    ggplot2::geom_point(ggplot2::aes(x = input_x, y = input_y),
                        alpha = 0.5,
                        size = 1.5) +
    # Add trend line with confidence interval
    ggplot2::geom_smooth(ggplot2::aes(x = input_x, y = input_y),
                         method = "lm",           # Linear regression
                         formula = y ~ x,         # Formula
                         se = TRUE,               # Show confidence interval
                         color = "royalblue",          # Trend line color
                         fill = "lightblue",      # Confidence interval fill color
                         alpha = 0.3,             # Transparency
                         linewidth = 0.8) +       # Line width

    # Labels
    ggplot2::labs(
      title = title_text,
      x = x_text,
      y = y_text
    ) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::theme(
      legend.position = "none",
      axis.title = ggplot2::element_text(size = text_size),
      axis.text = ggplot2::element_text(size = text_size * 0.9),
      axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5, vjust = 0.5),
      plot.title = ggplot2::element_text(size = text_size * 1.2, face = "bold",
                                         margin = ggplot2::margin(b = 3)),
      plot.title.position = "panel"
    )

  return(p)
}
