#!/bin/bash
set -e  # if error, stop

echo "=== START DATA TRANSFER ==="

# Move canopy height from Paul et al. 2026
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  --exclude='*_deg/' \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/canopy_height_paul_2026/ \
  /data_2/scratch/ting/veg_topo_data/data_raw/canopy_height_paul_2026/

# Move all processed data
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  --exclude='*_deg/' \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/ \
  /data_2/scratch/ting/veg_topo_data/data/

rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/global_rf_30m \
  /data_2/scratch/ting/veg_topo_data/data/

# Paul tree height
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data_raw/canopy_height_paul_2026 \
  /data_2/scratch/ting/veg_topo_data/data_raw/

# data/figures
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/homefs/tt22k003/veg_topo/data/figures/ \
  /data_2/scratch/ting/veg_topo_data/data/figures/

rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/homefs/tt22k003/veg_topo/data/figures/ \
  ~/veg_topo/data/figures/

# data/figures
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/homefs/tt22k003/veg_topo/data/figures/ \
  ~/veg_topo/data/figures/

# TWI 450m clean
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/global_twi_450m_clean \
  /data_2/scratch/ting/veg_topo_data/data/

# Veg 450m, exclude the sub-folder of splited tiles
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  --exclude='*/' \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/global_vegh_fveg_450m/ \
  /data_2/scratch/ting/veg_topo_data/data/global_vegh_fveg_450m/

# sw_in_450m, exclude the sub-folder of splited tiles
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  --exclude='*/' \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/global_sw_in_450m/ \
  /data_2/scratch/ting/veg_topo_data/data/global_sw_in_450m/

# global_dem_slope_aspect_450m, exclude the sub-folder of splited tiles
rsync -avhP --no-perms --no-owner --no-group \
  --exclude='*/' --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/global_dem_slope_aspect_450m/ \
  /data_2/scratch/ting/veg_topo_data/data/global_dem_slope_aspect_450m/

# global_r_H_R_5km
rsync -avhP --no-perms --no-owner --no-group \
  --exclude='*/' --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/global_r_H_R_5km/ \
  /data_2/scratch/ting/veg_topo_data/data/global_r_H_R_5km/

# global_cor_twi_vegh
rsync -avhP --no-perms --no-owner --no-group \
  --exclude='*/' --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/global_cor_twi_vegh/ \
  /data_2/scratch/ting/veg_topo_data/data/global_cor_twi_vegh/

# reg_sample
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/reg_sample/ \
  /data_2/scratch/ting/veg_topo_data/data/reg_sample/

# reg_validate_500m
rsync -avhP --no-perms --no-owner --no-group \
  --human-readable -i --info=progress2 \
  tt26y485@submit03.unibe.ch:/storage/scratch/giub_geco/tting/data/reg_validate_500m/ \
  /data_2/scratch/ting/veg_topo_data/data/reg_validate_500m/


echo "=== ALL TRANSFERS COMPLETED ==="
