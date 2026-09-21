#!/bin/bash
set -e  # if error, stop

echo "=== START DATA TRANSFER ==="
# All data in scratch folder.
rsync --human-readable -i --info=progress2 -avhP --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/ting/data/

rsync --human-readable -i --info=progress2 -avhP --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data_raw \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# ========= Raw data ===========
# VEGH 10m
rsync --human-readable -i --info=progress2 -avhP \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  --no-owner --no-group \
  /data_2/archive/vegheight_lang_2023 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# DEM 30m
rsync --human-readable -i --info=progress2 -avhP --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data_raw/COP30_OpenTopography \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# TWI 30m -- Ho et al. 2025
rsync -avhP --info=progress2 --no-owner --no-group \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  /data_2/archive/twi_ho_2025 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# TWI 450m -- Marthew et. al 2015
rsync -avhP --info=progress2 --no-owner --no-group \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  /data/archive/gti_marthews_2015 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# MI 950m -- Zomer 2022
rsync -avhP --info=progress2 --no-owner --no-group \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  /data/archive/aridityindex_zomer_2022 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# Climatic Variables
rsync -avhP --info=progress2 --no-owner --no-group \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  /data/archive/worldclim_fick_2017 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# Koeppen-Geiger
rsync --human-readable -i --info=progress2 -avhP \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  --no-owner --no-group \
  /data/archive/koeppengeiger_beck_2018 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# CCI land cover 300m
rsync --human-readable -i --info=progress2 -avhP \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  --no-owner --no-group \
  /data/archive/landcover_defourny_2023 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# WDPA: world database on Protected area
rsync --human-readable -i --info=progress2 -avhP --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data_raw/wdpa_2025 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# Ecoregions (biomes)
rsync --human-readable -i --info=progress2 -avhP --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data_raw/ecoregion2017 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# Airborne Lidar data using in Lang et al.(2023)
rsync --human-readable -i --info=progress2 -avhP --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data_raw/airborne_lidar_lang_2023 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/

# ========= Processed data ===========

# TWI (clean) 450m
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_twi_450m_clean \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# flc 5km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_flc_5km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# r_H_TWI 5km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_cor_twi_vegh \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# r_H_Rin 5km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_r_H_R_5km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# k-means 5km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_kmeans_5km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# mi 5km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_mi_5km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# mat 5km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_mat_5km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# VEGH 450m
rsync --human-readable -i --info=progress2 -avhP --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_vegh_fveg_450m \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# reg_sample
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  /data_2/scratch/ting/veg_topo_data/data/reg_sample \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# global_dem_slope_aspect_450m
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  --chmod=Du+rwx,Dgo+rx,Fu+rw,Fgo+r \
  /data_2/scratch/ting/veg_topo_data/data/global_dem_slope_aspect_450m \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# ecoregion_5km
rsync --human-readable -i --info=progress2 -avhP --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_ecoregion_5km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# dem 55km: dem_sd_55km, dem_rg_55km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_dem_55km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# flc 55km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_flc_55km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# mi 55km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_mi_55km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# fpa 55km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_fpa_55km_2025 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# fveg 55km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_fveg_55km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

# mat 55km
rsync -avhP --info=progress2 --no-perms --no-owner --no-group \
  /data_2/scratch/ting/veg_topo_data/data/global_mat_55km \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/

echo "=== ALL TRANSFERS COMPLETED ==="
