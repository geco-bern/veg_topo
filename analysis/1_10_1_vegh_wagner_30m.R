
# ============= Setup =======================================================
library(terra)
source(here::here("R/config.R"))

# Set workers based on system
if (hostname == "dash") workers <- 16 else workers <- 64
message("→ Using ", workers, " workers")

# Create temp directory
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

# Create output directory
output_dir <- dirname(vegh_wagner_30m_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("→ Output directory: ", output_dir)

# Configure temp paths for terra and GDAL
terra::terraOptions(tempdir = temp_dir)
Sys.setenv(TMPDIR = temp_dir, GDAL_TEMP_DIR = temp_dir, TEMP = temp_dir, TMP = temp_dir)
message("→ Temp directory: ", temp_dir)

tictoc::tic()
# ============= Prepare VRT using file list ==================================

# Get all tile files
# Found 22063 tile files (~185GB)
tile_files <- list.files(vegh_wagner_dir, pattern = "PAGE_.*_height_mean\\.tif$",
                         full.names = TRUE)[1:10] #  for test add [1:100]!!!!!
cat("→ Found", length(tile_files), "tile files (~185GB)\n")

# Create file list in output directory
list_file <- file.path(output_dir, "tile_list.txt")
writeLines(tile_files, list_file)
message("→ File list created: ", list_file)

# Create VRT with NoData handling using file list
message("→ Creating VRT...")
vrt_file <- file.path(vegh_wagner_dir, "all_tiles.vrt")

system2("gdalbuildvrt",
        args = c("-overwrite",
                 "-srcnodata", "0",
                 "-vrtnodata", "0",
                 "-input_file_list", shQuote(list_file),
                 shQuote(vrt_file)),
        stdout = TRUE, stderr = TRUE)

# ============= Warp (Reproject + Resample) ===================================

message("→ Starting warp (45-90 min)...")
cat("→ Using", workers, "cores\n")

# Remove output if it exists
if (file.exists(vegh_wagner_30m_path)) unlink(vegh_wagner_30m_path)

system2("gdalwarp",
        args = c(

          # Target CRS and resolution
          "-t_srs", "EPSG:4326",
          "-tr", "0.00025", "0.00025",

          # Resampling method
          "-r", "lanczos",

          # NoData handling
          "-srcnodata", "0",
          "-dstnodata", "0",

          # Memory configuration
          "-wm", "131072",
          "--config", "GDAL_CACHEMAX", "262144",
          "--config", "GDAL_TEMP_DIR", temp_dir,

          # Parallel processing
          "-multi",
          "-wo", paste0("NUM_THREADS=", workers),

          # Output compression
          "-co", "COMPRESS=ZSTD",
          "-co", "ZSTD_LEVEL=1",
          "-co", "BIGTIFF=YES",

          # Input/Output
          shQuote(vrt_file),
          shQuote(vegh_wagner_30m_path)
        ),
        stdout = TRUE, stderr = TRUE)

tictoc::toc()

# ============= WVerify =====================================================

if (file.exists(vegh_wagner_30m_path)) {
  file_size <- round(file.info(vegh_wagner_30m_path)$size / 1e9, 1)
  cat("→ Success!\n")
  cat("→ Output:", vegh_wagner_30m_path, "\n")
  cat("→ Size:", file_size, "GB\n")

  # Optional: clean up
  unlink(temp_dir, recursive = TRUE)
} else {
  cat("✗ Error: Output not created\n")
  # quit(status = 1)
}
