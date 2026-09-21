
# ----Set-Up----
library(terra)
library(tidyr)
library(dplyr)
library(car) # Box-Cox Transformation
library(ClusterR) # MiniBatchKmeans

library(ggplot2)
library(patchwork)
library(tidyterra)
library(sf)

source(here::here("R/config.R"))
source(here::here("R/plot_density_grid.R"))

# Create output directory
if (!dir.exists(dirname(kmeans_map_8c_path))) dir.create(dirname(kmeans_map_8c_path), recursive = TRUE)

# ----Data-Load-----
message("Data pre...")
# Load resampled raster datasets (AI, TWI, fused)
mi_5km_r <- terra::rast(mi_5km_file) * 0.0001 # Multiply all values by 0.0001 to get original value
fused_5km_r <- terra::rast(fused_5km_file)
cor_twi_vegh_5km_r <- terra::rast(cor_twi_vegh_mosaic_file)

# Stack rasters into a single SpatRaster
stacked <- c(cor_twi_vegh_5km_r,
             fused_5km_r,
             mi_5km_r)

# Convert to data frame for k-means clustering
df <- as.data.frame(stacked, xy = TRUE, na.rm = TRUE)
colnames(df) <- c("lon", "lat", "cor", "fused", "mi")

# save crs for write NetCDF document
tar_crs = terra::crs(mi_5km_r)
rm(mi_5km_r, fused_5km_r, stacked); gc()

# ----Data-Raw-Density----
# Select columns
cols <- c("mi", "cor", "fused")
plot_density_grid(
  df,
  cols,
  nrow = 1,
  width = 6,
  height = 2,
  save_path = here::here("data/figures/4_03_1_kmeans_data_raw_ds.png")
)

# ----Data-Pre----

# Use Box-Cox transformation to select the optimal lambda (power)
# This helps make the 'mi' variable more normally distributed
pt <- powerTransform(df$mi)
lambda <- pt$lambda
cat("Optimal lambda:", lambda, "\n")

# Apply the Box-Cox transformation using the selected lambda
df$mi_trans <- bcPower(df$mi, lambda)

# Standardize the input variables for k-means clustering
# This centers the data (mean = 0) and scales it (SD = 1)
df_k <- as.data.frame(scale(df[, c("cor", "fused", "mi_trans")]))

# ----Density-Data-Pre----

# Select columns
cols <- c("cor", "fused","mi_trans")
plot_density_grid(
  df_k,
  cols,
  main_title = "Standardized variables density distribution",
  nrow = 1,
  width = 7,
  height = 2.5,
  save_path = here::here("data/figures/4_03_1_kmeans_data_pre_ds.png")
)

# ----K-Means----

message("K-means...")
tictoc::tic()
set.seed(123)
# Use Mini-Batch K-Means
mb_km <- MiniBatchKmeans(
  df_k,                # scaled data
  clusters = 8,         # number of clusters
  batch_size = 10000,   # size of each mini-batch
  num_init = 10,        # number of random initializations
  max_iters = 100       # max iterations
)

# Add clustering results to data frame
df$cluster8c <- predict(mb_km, df_k)
tictoc::toc()

rm(mb_km, df_k); gc()

# ----Save-Cluster-Raster----

message("Save k-means result...")
# Reload an original raster as template
template_raster <- terra::rast(cor_twi_vegh_5km_r)

# Create empty raster with same properties as template
kmeans_8c_r <- terra::rast(template_raster)
values(kmeans_8c_r) <- NA  # Initialize with NA

# Convert clustering results to spatial points
df_sf <- sf::st_as_sf(df[, c("lon", "lat", "cluster8c")],
                      coords = c("lon", "lat"),
                      crs = tar_crs)

# Rasterize the points
kmeans_8c_r <- terra::rasterize(df_sf, template_raster, field = "cluster8c")

# Write to NetCDF
terra::writeCDF(kmeans_8c_r,
                filename = kmeans_map_8c_path,
                names = "cluster8c",
                overwrite = TRUE)

message(paste0("Cluster map saved to: ", kmeans_map_8c_path))

kmeans_map_8c_path
