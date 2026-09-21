
# ~ 128 min: UBELIX, 16 cores

# ------Set up------------------------------------------------------------------
# Load packages
library(terra)
library(furrr)
library(future)

source(here::here("R/config.R"))
source(here::here("R/raster_preprocess_save.R")) # Set Na value and aggregation

# Set worker numbers for different system
if (hostname == "dash") workers = 4 else workers = 16
message("→ using ", workers, " workers")

# -------Configuration----------------------------------------------------------

# Input: read all 10m vegetation data paths (2651)
vegh_10m_tiles_path <- fs::dir_ls(path = vegh_10m_tiles_dir, glob = "*_Map.tif")
vegh_10m_tiles_path_sub <- vegh_10m_tiles_path[1:2651]  # subset for testing (limit number of tiles)!!!!!!

# Output: 2651 tiles with 450m resolution
if (!dir.exists(vegh_450m_tiles_dir)) {
  dir.create(vegh_450m_tiles_dir, recursive = TRUE)
}

# target layer
res_tar <- c(15/3600, 15/3600)

# -----Aggregation: 10m -> 450m ------------------------------------------------

agg_vegh_fveg <- function(vegh_10m_file) {

  tryCatch({

    # --- tile info ---
    tile_id <- basename(vegh_10m_file)
    tictoc::tic(paste0("Processing tile: ", tile_id))
    t0 <- Sys.time()

    # set output file
    vegh_450m_file <- file.path(
      vegh_450m_tiles_dir,
      paste0(sub("\\.tif$", "", basename(vegh_10m_file)), "_to450m.nc")
    )

    # set output file
    fveg_450m_file <- file.path(
      vegh_450m_tiles_dir,
      paste0(sub("\\.tif$", "", basename(vegh_10m_file)), "_to450m_fveg.nc")
    )

    # Check if files have been processed
    if (fs::file_exists(vegh_450m_file) && fs::file_exists(fveg_450m_file)) {
      return(TRUE)
    }

    # --- aggregate vegetation height ---
    # Set 0 as NA value (0m canopy height represents not vegetated or water according to Lang et al. (2019))
    # Aggregates using TWI data from Marthews et al. (2015) as target resolution
    # Return the saved path
    agg_vegh <- suppressMessages(
      raster_preprocess_save(
      input = vegh_10m_file,
      output = vegh_450m_file,
      res_tar = res_tar,
      na_value = 0,
      fun = mean,
      varname = "vegh",
      if_aggregate = TRUE,
      if_round_fact = TRUE,
      if_resample = FALSE,
      if_mask = FALSE,
      if_return_raster = FALSE
      ))

    # --- Aggregate fraction of vegetated area ---
    # Compute the fraction of pixels with vegetation (height > 0)
    # while treating 0 and NA as non-vegetated.
    # Aggregation is done at the target resolution defined by TWI data (Marthews et al., 2015).
    agg_fveg <- suppressMessages(
      raster_preprocess_save(
        input        = vegh_10m_file,
        output       = fveg_450m_file,
        res_tar      = res_tar,
        varname      = "fveg",
        if_aggregate = TRUE,
        fun = function(x, na.rm) {
          total <- length(x)                  # total number of pixels INCLUDING NA
          if (all(is.na(x))) return(NA)       # if all NA return NA
          veg_count <- sum(!is.na(x) & x > 0)  # count of vegetated pixels
          return(veg_count / total)            # fraction over TOTAL pixels, including NA
        },
        if_round_fact = TRUE,
        if_resample = FALSE,
        if_mask = FALSE,
        if_return_raster = FALSE
      ))

    # --- Print proccess time ---
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("Tile %s completed [%.1f mins]", tile_id, elapsed_mins))
    tictoc::toc()

    return(TRUE)

  }, error = function(e) {
    elapsed_mins <- difftime(Sys.time(), t0, units = "mins")
    message(sprintf("❌ Tile %s failed after %.1f mins: %s", basename(vegh_10m_file), elapsed_mins, e$message))
    return(FALSE)
  })
}

# ----------------- Parallel execution for all tiles-----------------
# Set up cluster plan
plan(cluster, workers = workers)

# Run in parallel
tictoc::tic("🚀 Parallel processing of tiles")

results <- future_map(
  vegh_10m_tiles_path_sub,
  agg_vegh_fveg,
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
# r0 <- rast(file.path(vegh_10m_tiles_dir,"ETH_GlobalCanopyHeight_10m_2020_N00E033_Map.tif"))
# r1 <- rast("/storage/scratch/giub_geco/tting/data/global_vegh_fveg_450m/3_3_deg/ETH_GlobalCanopyHeight_10m_2020_N00E033_Map_to450m.nc")
# r2 <- rast("/data_2/scratch/ting/veg_topo_data/data/global_vegh_450m/..._to450m_fveg.nc")
# plot(r0, xlim=c(6.6,6.7), ylim=c(0.1,0.2))
# plot(r1, xlim=c(6.6,6.7), ylim=c(0.1,0.2))
# plot(r2, xlim=c(6.6,6.7), ylim=c(0.1,0.2))
