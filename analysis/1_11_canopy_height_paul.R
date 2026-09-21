# this script just for check the downloaded files
library(terra)

source(here::here("R/config.R"))
file <- "/storage/scratch/giub_geco/tting/data_raw/canopy_height_paul_2026/ECHOSAT_T01G.tif"
r <- rast(file)
r

file <- "/storage/scratch/giub_geco/tting/data_raw/canopy_height_paul_2026/ECHOSAT_T01G_2020.tif"
r <- rast(file)
r

file <- "/storage/scratch/giub_geco/tting/data_raw/canopy_height_paul_2026/ECHOSAT_T01U.tif"
r <- rast(file)
r
plot(r)


file <- "/storage/scratch/giub_geco/tting/data_raw/canopy_height_paul_2026/ECHOSAT_T05U_2020.tif"
r <- rast(file)
r
summary(r)
plot(r)
