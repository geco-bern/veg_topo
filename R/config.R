# ----------------------- NOTES ------------------------
# - ext_global requires the terra package to work!!!
#
# - project_root:
#   Source code and figure output directory.
#
# - veg_topo_extr_dir:
#   Storage for large datasets (raw, intermediate, outputs).
#
# - Raw data locations:
#   - /data/archive/
#   - /data_2/archive/
#   - veg_topo_extr_dir/data_raw
#
# - Intermediate / output data:
#   - veg_topo_extr_dir/data/global_*   # global-scale products
#   - veg_topo_extr_dir/data/reg_*      # regional-scale products

# ----------------------- Raw data ------------------------
# Global extent (terra::ext)
ext_global <- tryCatch(
  ext(-180, 180, -60, 85),
  error = function(e) {
    message("⚠️ Failed to create extent: ", e$message)
    NULL
  }
)

# Auto-detect environment and set paths
hostname <- trimws(tolower(system("hostname", intern = TRUE)))

if (hostname == "dash") {
  # --- Workstation2: Base directories ---
  message("💻 Workstation detected: dash")
  project_root <- "~/veg_topo"
  veg_topo_extr_dir <- file.path("/data_2/scratch/ting/veg_topo_data")
  temp_dir <- file.path("/data_2/scratch/ting/veg_topo_data/tmp")

  # --- Workstation2: Raw data paths ---
  # Topographic Wetness Index (TWI)
  twi_450m_path <- file.path("/data/archive/gti_marthews_2015/data/ga2.nc")   # ~450 m
  twi_30m_path  <- file.path("/data_2/archive/twi_ho_2025/data/twi_edtm_m_30m_v20241230.tif")  # 30 m

  # Vegetation canopy height (10 m tiles)
  vegh_10m_tiles_dir <- file.path("/data_2/archive/vegheight_lang_2023/data/3deg_cogs")

  # Elevation (Copernicus DEM, ~30 m)
  dem_30m_copernicus_dir <- file.path(veg_topo_extr_dir, "data_raw/copernicus_dem_30m/copernicus_dem_30m")
  COP30_dir <- file.path(veg_topo_extr_dir, "data_raw/COP30_OpenTopography/COP30_hh")

  # Aridity / Moisture Index
  mi_950m_file <- file.path("/data/archive/aridityindex_zomer_2022/data/Global-AI_ET0_v3_annual/ai_v3_yr.tif")

  # WorldClim 2
  worldclim_1km_dir <- file.path("/data/archive/worldclim_fick_2017/data/")

  # ESA CCI Land Cover (300 m, 2020)
  cci_landcover_path <- file.path("/data/archive/landcover_defourny_2023/data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc")

  # Ecoregions (biomes)
  ecoregion_path <- file.path(veg_topo_extr_dir, "data_raw/ecoregion2017/Ecoregions2017/Ecoregions2017.shp")

  # Protected areas (WDPA, split shapefiles)
  pa_shp0 <- file.path(veg_topo_extr_dir, "data_raw/wdpa_2025/WDPA_WDOECM_Jul2025_Public_all_shp/WDPA_WDOECM_Jul2025_Public_all_shp_0/WDPA_WDOECM_Jul2025_Public_all_shp-polygons.shp")
  pa_shp1 <- file.path(veg_topo_extr_dir, "data_raw/wdpa_2025/WDPA_WDOECM_Jul2025_Public_all_shp/WDPA_WDOECM_Jul2025_Public_all_shp_1/WDPA_WDOECM_Jul2025_Public_all_shp-polygons.shp")
  pa_shp2 <- file.path(veg_topo_extr_dir, "data_raw/wdpa_2025/WDPA_WDOECM_Jul2025_Public_all_shp/WDPA_WDOECM_Jul2025_Public_all_shp_2/WDPA_WDOECM_Jul2025_Public_all_shp-polygons.shp")

  # Köppen–Geiger climate classification
  kg_present_0p0083_file <- file.path("/data/archive/koeppengeiger_beck_2018/data/Beck_KG_V1_present_0p0083.tif")  # ~1 km
  kg_present_0p083_file  <- file.path("/data/archive/koeppengeiger_beck_2018/data/Beck_KG_V1_present_0p083.tif")   # ~10 km
  kg_legend_file         <- file.path("/data/archive/koeppengeiger_beck_2018/data/legend.txt")

  # Airborne Lidar ALS and LVIS
  lidar_asl_dir <- file.path(veg_topo_extr_dir, "data_raw/airborne_lidar_lang_2023/airborne_lidar/ALS_MaxGEDIFootprint_GSD10m")
  lidar_lvis_dir <- file.path(veg_topo_extr_dir, "data_raw/airborne_lidar_lang_2023/airborne_lidar/LVIS_RH98_GSD10m")

  # VegH Wagner 2024
  vegh_wagner_dir <- file.path(veg_topo_extr_dir, "data_raw/canopy_height_wagner_2024/v1")
} else {

  # --- UBELIX: Base directories ---
  message("🖥️ HPC environment detected: ", hostname)
  project_root <- "~/veg_topo"
  veg_topo_extr_dir <- file.path("/storage/scratch/giub_geco/tting")
  temp_dir <- file.path("/storage/scratch/giub_geco/tting/tmp")

  # --- UBELIX: Raw data paths ---
  # Topographic Wetness Index (TWI)
  twi_450m_path <- file.path(veg_topo_extr_dir, "data_raw/gti_marthews_2015/data/ga2.nc")   # ~450 m
  twi_30m_path  <- file.path(veg_topo_extr_dir, "data_raw/twi_ho_2025/data/twi_edtm_m_30m_v20241230.tif")  # 30 m

  # Vegetation canopy height (10 m tiles)
  vegh_10m_tiles_dir <- file.path(veg_topo_extr_dir, "data_raw/vegheight_lang_2023/data/3deg_cogs")

  # Elevation (Copernicus DEM, ~30 m)
  dem_30m_copernicus_dir <- file.path(veg_topo_extr_dir, "data_raw/copernicus_dem_30m/copernicus_dem_30m")
  COP30_dir <- file.path(veg_topo_extr_dir, "data_raw/COP30_OpenTopography/COP30_hh")

  # Aridity / Moisture Index
  mi_950m_file <- file.path(veg_topo_extr_dir, "data_raw/aridityindex_zomer_2022/data/Global-AI_ET0_v3_annual/ai_v3_yr.tif")

  # WorldClim 2
  worldclim_1km_dir <- file.path(veg_topo_extr_dir, "data_raw/worldclim_fick_2017/data/")

  # ESA CCI Land Cover (300 m, 2020)
  cci_landcover_path <- file.path(veg_topo_extr_dir, "data_raw/landcover_defourny_2023/data/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.nc")

  # Ecoregions (biomes)
  ecoregion_path <- file.path(veg_topo_extr_dir, "data_raw/ecoregion2017/Ecoregions2017/Ecoregions2017.shp")

  # Protected areas (WDPA, split shapefiles)
  pa_shp0 <- file.path(veg_topo_extr_dir, "data_raw/wdpa_2025/WDPA_WDOECM_Jul2025_Public_all_shp/WDPA_WDOECM_Jul2025_Public_all_shp_0/WDPA_WDOECM_Jul2025_Public_all_shp-polygons.shp")
  pa_shp1 <- file.path(veg_topo_extr_dir, "data_raw/wdpa_2025/WDPA_WDOECM_Jul2025_Public_all_shp/WDPA_WDOECM_Jul2025_Public_all_shp_1/WDPA_WDOECM_Jul2025_Public_all_shp-polygons.shp")
  pa_shp2 <- file.path(veg_topo_extr_dir, "data_raw/wdpa_2025/WDPA_WDOECM_Jul2025_Public_all_shp/WDPA_WDOECM_Jul2025_Public_all_shp_2/WDPA_WDOECM_Jul2025_Public_all_shp-polygons.shp")

  # Köppen–Geiger climate classification
  kg_present_0p0083_file <- file.path(veg_topo_extr_dir, "data_raw/koeppengeiger_beck_2018/data/Beck_KG_V1_present_0p0083.tif")  # ~1 km
  kg_present_0p083_file  <- file.path(veg_topo_extr_dir, "data_raw/koeppengeiger_beck_2018/data/Beck_KG_V1_present_0p083.tif")   # ~10 km
  kg_legend_file         <- file.path(veg_topo_extr_dir, "data_raw/koeppengeiger_beck_2018/data/legend.txt")

  # Airborne Lidar ALS and LVIS
  lidar_asl_dir <- file.path(veg_topo_extr_dir, "data_raw/airborne_lidar_lang_2023/airborne_lidar/ALS_MaxGEDIFootprint_GSD10m")
  lidar_lvis_dir <- file.path(veg_topo_extr_dir, "data_raw/airborne_lidar_lang_2023/airborne_lidar/LVIS_RH98_GSD10m")
}

# ----------------------- Processed data ------------------------

# 1_01 Topographic Wetness Index (cleaned)
twi_450m_mosaic_clean_path <- file.path(veg_topo_extr_dir,"data/global_twi_450m_clean/ga2_clean.tif")

# 1_02 Vegetation height and fraction of vegetated area
vegh_450m_tiles_dir   <- file.path(veg_topo_extr_dir,"data/global_vegh_fveg_450m/3_3_deg")
vegh_450m_mosaic_path <- file.path(veg_topo_extr_dir,"data/global_vegh_fveg_450m/vegh_450m_2020_mosaic.tif")
fveg_real_450m_mosaic_path <- file.path(veg_topo_extr_dir,"data/global_vegh_fveg_450m/fveg_real_450m_2020_mosaic.tif")
fveg_real_55km_path   <- file.path(veg_topo_extr_dir,"data/global_fveg_55km/fveg_real_55km.nc")
fveg_55km_path        <- file.path(veg_topo_extr_dir,"data/global_fveg_55km/fveg_55km.nc")

# 1_03 Elevation (DEM, slope, aspect at ~450 m)
dem_450m_tiles_dir    <- file.path(veg_topo_extr_dir,"data/global_dem_slope_aspect_450m/dem_1_1_deg")
slope_450m_tiles_dir  <- file.path(veg_topo_extr_dir,"data/global_dem_slope_aspect_450m/slope_1_1_deg")
aspect_450m_tiles_dir <- file.path(veg_topo_extr_dir,"data/global_dem_slope_aspect_450m/aspect_1_1_deg")

dem_450m_mosaic_path    <- file.path(veg_topo_extr_dir,"data/global_dem_slope_aspect_450m/dem_450m.tif")
slope_450m_mosaic_path  <- file.path(veg_topo_extr_dir,"data/global_dem_slope_aspect_450m/slope_450m.tif")
aspect_450m_mosaic_path <- file.path(veg_topo_extr_dir,"data/global_dem_slope_aspect_450m/aspect_450m.tif")

# Shortwave radiation (terrain vs flat)
# sw_in_uneven_450m_tile_dir              <- file.path(veg_topo_extr_dir,"data/global_sw_in_450m/uneven_1_1_deg")
# sw_in_uneven_450m_path           <- file.path(veg_topo_extr_dir,"data/global_sw_in_450m/sw_in_uneven_450m.tif")   # terrain
# sw_in_flat_450m_tile_dir              <- file.path(veg_topo_extr_dir,"data/global_sw_in_450m/flat_1_1_deg")
# sw_in_flat_450m_path             <- file.path(veg_topo_extr_dir,"data/global_sw_in_450m/sw_in_flat_450m.tif")     # flat
# sw_in_terrain_effect_450m_path   <- file.path(veg_topo_extr_dir,"data/global_sw_in_450m/sw_in_terrain_effect_450m.tif")  # ratio
rin_450m_tiles_dir <- file.path(veg_topo_extr_dir,"data/global_rin_450m/1_1_deg")  # ratio
rin_450m_mosaic_path   <- file.path(veg_topo_extr_dir,"data/global_rin_450m/radiation_index_450m.tif")  # ratio

# Elevation variability (~55 km)
dem_sd_55km_path            <- file.path(veg_topo_extr_dir,"data/global_dem_55km/dem_sd_55km.nc")
dem_rg_95p_05p_55km_path    <- file.path(veg_topo_extr_dir,"data/global_dem_55km/dem_rg_95p_05p_55km.nc")  # 95th–05p percentile

# 1_04 Moisture index
mi_5km_file  <- file.path(veg_topo_extr_dir,"data/global_mi_5km/mi_5km.nc")
mi_55km_file <- file.path(veg_topo_extr_dir,"data/global_mi_55km/mi_55km.nc")

# 1_05 Land cover fractions
flc_tile_dir   <- file.path(veg_topo_extr_dir,"data/global_flc_5km/30_30_deg")
fused_5km_file <- file.path(veg_topo_extr_dir,"data/global_flc_5km/fused_5km.nc")
fbare_5km_file <- file.path(veg_topo_extr_dir,"data/global_flc_5km/fbare_5km.nc")
fwater_5km_file<- file.path(veg_topo_extr_dir,"data/global_flc_5km/fwater_5km.nc")
fsnow_5km_file <- file.path(veg_topo_extr_dir,"data/global_flc_5km/fsnow_5km.nc")

fused_55km_file <- file.path(veg_topo_extr_dir,"data/global_flc_55km/fused_55km.nc")
fbare_55km_file <- file.path(veg_topo_extr_dir,"data/global_flc_55km/fbare_55km.nc")
fwater_55km_file<- file.path(veg_topo_extr_dir,"data/global_flc_55km/fwater_55km.nc")
fsnow_55km_file <- file.path(veg_topo_extr_dir,"data/global_flc_55km/fsnow_55km.nc")

# 1_06 Climatic variables (1970–2000)
mat_450m_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_mat_450m_1970_2000.tif")
mat_1km_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_mat_1km_1970_2000.tif")
mat_5km_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_mat_5km_1970_2000.nc")
mat_55km_file <- file.path(veg_topo_extr_dir,"data/global_climvar_55km/global_mat_55km_1970_2000.nc")

map_450m_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_map_450m_1970_2000.tif")
map_1km_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_map_1km_1970_2000.tif")
map_5km_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_map_5km_1970_2000.nc")
map_55km_file <- file.path(veg_topo_extr_dir,"data/global_climvar_55km/global_map_55km_1970_2000.nc")

srad_450m_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_srad_450m_1970_2000.tif")
srad_1km_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_srad_1km_1970_2000.tif")
srad_5km_file  <- file.path(veg_topo_extr_dir,"data/global_climvar_5km/global_srad_5km_1970_2000.nc")
srad_55km_file <- file.path(veg_topo_extr_dir,"data/global_climvar_55km/global_srad_55km_1970_2000.nc")

# 1_07 Biomes (ecoregions)
ecoregion_5km_path <- file.path(veg_topo_extr_dir,"data/global_ecoregion_5km/Ecoregions2017_5km.nc")

# 1_08 Protected areas fraction
fpa_55km_path <- file.path(veg_topo_extr_dir,"data/global_fpa_55km_2025/WDPA_WDOECM_Jul2025_55km.nc")

# 1_10 Amazon vegetation height data (Wagner et. al. 2024)
vegh_wagner_30m_path <- file.path(veg_topo_extr_dir,"data/vegh_wagner_30m/vegh_2015-12-01_END_2024-06-01_height_mean_wagner_30m.tif")
vegh_wagner_450m_path <- file.path(veg_topo_extr_dir,"data/vegh_wagner_450m/vegh_2015-12-01_END_2024-06-01_height_mean_wagner_450m.tif")
# ----------------------- Analysis products ------------------------

# 2_01 Valid tiles (global analysis)
valid_tiles_info_path <- file.path(project_root,"data/valid_tiles_info.rds")
valid_tiles_vect_path <- file.path(project_root,"data/valid_tiles.gpkg")

# 2_02 Correlation (H ~ TWI)
cor_twi_vegh_tiles_dir       <- file.path(veg_topo_extr_dir,"data/global_cor_twi_vegh/30_30_deg")
cor_twi_vegh_mosaic_file     <- file.path(dirname(cor_twi_vegh_tiles_dir),"cor_twi_vegh_5km_mosaic.nc")
pval_cor_twi_vegh_mosaic_file<- file.path(dirname(cor_twi_vegh_tiles_dir),"pval_cor_twi_vegh_5km_mosaic.nc")

# 3_01 Correlation (H ~ R)
r_H_R_tiles_dir  <- file.path(veg_topo_extr_dir,"data/global_r_H_R_5km/30_30_deg")
r_H_R_5km_path   <- file.path(veg_topo_extr_dir,"data/global_r_H_R_5km/r_H_R_5km.nc")
pval_r_H_R_5km_path <- file.path(veg_topo_extr_dir,"data/global_r_H_R_5km/pval_r_H_R_5km.nc")

# 4_01 K-means clustering (8 classes)
kmeans_map_8c_path <- file.path(veg_topo_extr_dir,"data/global_kmeans_5km/kmeans_corth_fused_mi_5km_8c.nc")

# 5_01 Sampled regions
reg_sample_dir       <- file.path(veg_topo_extr_dir,"data/reg_sample")
reg_sample_info_path <- file.path(reg_sample_dir,"reg_sample_info.rds")
reg_sample_vect_path <- file.path(reg_sample_dir,"reg_sample_vect.gpkg")

# 5_02 Validation results (regional)
reg_validate_dir <- file.path(veg_topo_extr_dir,"data/reg_validate_500m")

# 6_01 Validation of Height data
h_validation_dir <- file.path(veg_topo_extr_dir,"data/h_validate_500m")

# 7_01 Sample data for RF
rf_sample_data_tiles_dir <- file.path(veg_topo_extr_dir,"data/global_rf_30m/1_1_deg_sample_data")
rf_models_dir <- file.path(veg_topo_extr_dir,"data/global_rf_30m/models")
valid_win_path <- file.path(veg_topo_extr_dir, "data/global_rf_30m/global_valid_win_5km.tif")
rf_tiles_path <- file.path(veg_topo_extr_dir, "data/global_rf_30m/global_rf_tiles.parquet")
