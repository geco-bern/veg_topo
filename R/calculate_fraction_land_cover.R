# ==============================================================================
# Function: calculate_fraction_land_cover
#
# Purpose:
#   Calculate fractions of different land-use types (fused, fbare, fwater, fsnow)
#   from high-resolution land cover data (5km) and optionally save as NetCDF.
#
# Inputs:
#   - df_win      : Data frame containing columns 'lon_mid', 'lat_mid', 'lccs_class'
#   - output_file : Optional path to save the resulting raster (NetCDF format)
#
# Output:
#   - Returns a data frame with aggregated fractions per grid cell
#   - All grid cells are preserved, missing classes become 0
#   - Optionally saves as NetCDF raster
# ==============================================================================

calculate_fraction_land_cover <- function(df_win, output_file = NULL) {

  # ---- Reference land cover classes ----
  ref_classes <- c(10, 11, 12, 20, 30, 40, 190, 200, 201, 202, 210, 220)

  # Get all unique grid cells first
  all_grids <- df_win |>
    distinct(lon_mid, lat_mid)

  # ---- Compute fractions ----
  df_flc <- df_win |>
    # ensure factor levels
    mutate(lccs_class = factor(lccs_class, levels = ref_classes)) |>
    # count occurrences per class in each grid cell
    count(lon_mid, lat_mid, lccs_class, .drop = TRUE) |>
    # compute proportion within each grid cell
    group_by(lon_mid, lat_mid) |>
    mutate(prop = n / sum(n)) |>
    # aggregate fractions for fused, fbare, fwater, fsnow
    summarise(
      fused  = sum(prop[lccs_class %in% c(10, 11, 12, 20, 190)], na.rm = TRUE) +
        0.75 * sum(prop[lccs_class == 30], na.rm = TRUE) +
        0.25 * sum(prop[lccs_class == 40], na.rm = TRUE),
      fbare  = sum(prop[lccs_class %in% c(200, 201, 202)], na.rm = TRUE),
      fwater = sum(prop[lccs_class %in% c(210)], na.rm = TRUE),
      fsnow  = sum(prop[lccs_class %in% c(220)], na.rm = TRUE),
      .groups = "drop"
    ) |>
    # Replace NA with 0 to preserve all grid cells
    mutate(across(c(fused, fbare, fwater, fsnow), ~ replace_na(.x, 0)))

  # Join back to ensure all grids are preserved
  df_flc <- all_grids |>
    left_join(df_flc, by = c("lon_mid", "lat_mid")) |>
    mutate(across(c(fused, fbare, fwater, fsnow), ~ replace_na(.x, 0)))

  # ---- Save as raster if output_file provided ----
  if (!is.null(output_file)) {

    flc_r <- terra::rast(
      df_flc[, c("lon_mid", "lat_mid", "fused", "fbare", "fwater", "fsnow")],
      type = "xyz",
      crs  = "EPSG:4326"
    )

    names(flc_r) <- c("fused", "fbare", "fwater", "fsnow")

    terra::writeCDF(flc_r, output_file, overwrite = TRUE)
  }

  return(df_flc)
}
