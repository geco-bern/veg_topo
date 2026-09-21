#' Plot Boxplot for Clustered Data
#'
#' This function generates a boxplot for a specified numeric variable
#' grouped by a categorical variable (cluster). Outliers are hidden
#' but y-axis range can be manually controlled.
#'
#' @param data Data frame containing the data.
#' @param xvar Character. Name of the categorical variable (cluster) for x-axis.
#' @param yvar Character. Name of the numeric variable to plot on y-axis.
#' @param ylab Character. Y-axis label. Default is NULL.
#' @param fill_colors Vector of colors for each cluster. Default is NULL.
#' @param show_legend Logical. If TRUE, show legend; otherwise hide legend. Default is FALSE.
#' @param text_size Numeric. Base text size. Default is 7.
#' @param ylim Numeric vector of length 2. Manual y-axis limits c(lower, upper).
#'        Default is NULL (automatically determined from data).
#'
#' @return A ggplot object of the boxplot.
#'
#' @examples
#' # Auto y-axis range
#' p <- plot_boxplot(data = df, xvar = "cluster8c", yvar = "mi", ylab = "MI")
#'
#' # Manual y-axis range
#' p <- plot_boxplot(data = df, xvar = "cluster8c", yvar = "mi", ylab = "MI",
#'                   ylim = c(-0.5, 0.5))
#' print(p)
plot_boxplot <- function(data, xvar, yvar, ylab = NULL, show_legend = FALSE,
                         text_size = 7, ylim = NULL) {

  # ---- 只去除 NA，保留所有数据（包括离群点） ----
  data_clean <- data[!is.na(data[[yvar]]), ]

  message(sprintf("Using %d rows (NA removed, outliers kept)", nrow(data_clean)))

  # ---- fill colors ----
  fill_colors <- setNames(
    c("#E78AC3", "#FC8D62", "#FFD92F", "#E5C494", "#B3B3B3", "#66C2A5", "#8DA0CB", "#A6D854"),
    cluster_labels)

  # ---- 创建基础图 ----
  p <- ggplot(data_clean, aes(x = .data[[xvar]], y = .data[[yvar]], fill = .data[[xvar]])) +
    stat_boxplot(geom = "errorbar", width = 0.6, linewidth = 0.2) +
    geom_boxplot(
      width = 0.8,
      linewidth = 0.2,
      outlier.shape = NA,  # 不显示离群点
      coef = 1.5           # 离群点判定标准：IQR的1.5倍
    ) +
    { if (yvar == "cor") geom_hline(yintercept = 0, linetype = "dashed", color = "red") } +
    scale_fill_manual(values = fill_colors, name = "Group", guide = if (show_legend) "legend" else "none") +
    labs(title = NULL, y = ylab, x = "") +
    scale_x_discrete(drop = TRUE, expand = c(0.1, 0.1)) +
    theme_bw(base_size = text_size) +
    theme(
      legend.position = "none",
      legend.text = element_text(size = text_size),
      legend.title = element_text(size = text_size),
      axis.title = element_text(size = text_size),
      axis.text.x = element_blank(),
      axis.text = element_text(size = text_size),
      plot.title = element_blank()
    )

  # ---- 添加 y 轴范围控制 ----
  if (!is.null(ylim)) {
    # 检查 ylim 是否有效
    if (!is.numeric(ylim) || length(ylim) != 2) {
      stop("ylim must be a numeric vector of length 2: c(lower, upper)")
    }
    if (ylim[1] >= ylim[2]) {
      stop("ylim[1] must be less than ylim[2]")
    }

    # 使用 coord_cartesian 限制显示范围（不改变数据）
    p <- p + coord_cartesian(ylim = ylim)
  }

  return(p)
}
