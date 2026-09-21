# --------------- Setup --------------------------------------------------------
library(terra)
library(ggplot2)

source(here::here("R/config.R"))
source(here::here("R/mosaic_tiles.R"))
source(here::here("R/create_aligned_template.R"))
source(here::here("R/raster_preprocess_save.R"))


r <- rast("/data_2/scratch/ting/veg_topo_data/data_raw/airborne_lidar_lang_2023/airborne_lidar/ALS_MaxGEDIFootprint_GSD10m/32TMT.tif")
# r <- rast("/data_2/scratch/ting/veg_topo_data/data_raw/airborne_lidar_lang_2023/airborne_lidar/ALS_MaxGEDIFootprint_GSD10m/32TMS.tif")

# 提取所有值（向量形式）
vals <- values(r, matlab = FALSE)

# 计算?%分位数阈值
thresh <- quantile(vals, 0.99, na.rm = TRUE)

# 提取最高?%的数据
top <- vals[vals >= thresh]
summary(top)

df_top <- data.frame(
  x = 1:length(top),
  y = top
)

# 绘制hex图
ggplot(df_top, aes(x = x, y = y)) +
  geom_hex(bins = 100) +  # bins 控制六边形数量
  scale_fill_viridis_c() +  # 使用viridis颜色方案
  labs(title = "32TMT - Top ?% Canopy Height Hexbin Plot",
       subtitle = paste("?th percentile threshold:", round(thresh, 2), "m"),
       x = "Data Point Index (sorted by height)",
       y = "Height (m)",
       fill = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
