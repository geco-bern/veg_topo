#' Transform raster to target CRS/resolution: zonal (exact), aggregate + project (recommended for CRS change), or resample (same CRS).
#'
#' # Pipeline:
#   if_zonal      -> exact_zonal -> mask
#   if_project    -> (aggregate?) -> project -> mask
#   if_aggregate  -> aggregate -> (resample?) -> mask
#   if_resample   -> resample -> mask
#'
#' @param input Raster file or SpatRaster
#' @param output Output path (folder or file vector)
#' @param res_tar Target resolution c(x,y) for aggregate
#' @param target Template raster for zonal/project/resample/mask
#' @param varname NetCDF variable name (default "band")
#' @param na_value Value to set as NA
#' @param if_zonal Exact zonal aggregation to target grid (requires target)
#' @param if_aggregate Traditional aggregation to res_tar (default TRUE)
#' @param if_project Project to target CRS (implies resample)
#' @param if_resample Resample to target grid (same CRS only)
#' @param if_mask Mask to target extent
#' @param if_round_fact Round aggregation factors to integer (default TRUE)
#' @param if_return_raster Return SpatRaster (default TRUE)
#' @param fun Aggregation function (default mean)
#' @param overwrite Overwrite output files (default TRUE)
#'
#' @return SpatRaster or output file paths
#'
#' @examples
#' # Exact zonal aggregation (slowest but most accurate)
#' raster_preprocess("input.tif", target = "grid.tif", if_zonal = TRUE)
#'
#' # Aggregate then project (recommended for CRS change)
#' raster_preprocess("input.tif", target = "grid.tif", res_tar = 30,
#'                   if_aggregate = TRUE, if_project = TRUE)
#'
#' # Simple resample (same CRS)
#' raster_preprocess("input.tif", target = "grid.tif", if_resample = TRUE)
#'
#' @export
source(here::here("R/save_raster.R"))
raster_preprocess_save <- function(input,
                                   output = NULL,
                                   res_tar = NULL,
                                   target = NULL,
                                   varname = "band",
                                   na_value = NULL,
                                   if_project = FALSE,
                                   if_zonal = FALSE,
                                   if_aggregate = TRUE,
                                   fun = mean,
                                   if_resample = FALSE,
                                   if_round_fact = TRUE,
                                   if_mask = FALSE,
                                   if_return_raster = TRUE,
                                   overwrite = TRUE) {

  # --- Load input raster ---
  r_in <- if (is.character(input)) terra::rast(input) else input

  # Replace specified NA value with proper NA
  if (!is.null(na_value)) r_in[r_in == na_value] <- NA

  # --- Get target raster and resolution ---
  if (!is.null(target)) r_tar <- if (is.character(target)) terra::rast(target)[[1]] else target[[1]]
  if (!is.null(target) && is.null(res_tar)) res_tar <- c(terra::xres(r_tar), terra::yres(r_tar))
  if (is.null(res_tar))  warning("Both target and res_tar are NULL!")

  # ========== Aggregation method selection (mutually exclusive, zonal priority) ==========
  # --- Parameter validation ---
  if (if_zonal && if_aggregate) {
    warning("Both if_zonal and if_aggregate are TRUE. ",
            "Using if_zonal (priority), ignoring if_aggregate.")
  }

  if (if_zonal && is.null(target)) {
    stop("if_zonal = TRUE requires a 'target' raster.")
  }

  # --- Aggregation (Zonal/ Aggregate/ Disaggregate)---
  if (if_zonal && !is.null(target)) {

    # --- Method 1: Zonal Aggregate (precise, no interpolation, no alignment needed) ---

    message("Using ZONAL aggregate (high-res -> target grid)")

    if (if_aggregate) {
      message("  -> Note: if_aggregate = TRUE ignored because if_zonal = TRUE takes priority")
    }

    # Load or extract target raster
    r_tar <- if (is.character(target)) terra::rast(target)[[1]] else target[[1]]

    # Expand the input raster's spatial extent to align with the target grid
    # This ensures all target polygons have corresponding raster data
    # Areas outside original extent are set to NA
    # super important for fveg calculation!!!!
    r_in <- terra::extend(r_in, ext(r_tar))

    # Convert target raster to polygons (memory-efficient with chunked processing)
    target_polygons <- sf::st_as_sf(terra::as.polygons(r_tar,
                                                       dissolve = FALSE,
                                                       trunc = FALSE,
                                                       values = TRUE))

    # Efficient zonal extraction using exactextractr
    result <- exactextractr::exact_extract(
      r_in,
      target_polygons,
      fun = fun,
      progress = FALSE
    )

    # Create output raster with target's structure and new values
    r_out <- rast(r_tar)  # Create new raster from target's template
    terra::values(r_out) <- result  # Assign aggregated values

  } else if (if_aggregate && !is.null(res_tar)) {

    # --- Method 2: Traditional Aggregate (requires grid alignment) ---
    message("Using TRADITIONAL aggregate")

    # Calculate aggregation factors (input resolution -> target resolution)
    fact_x <- res_tar[1] / res(r_in)[1]
    fact_y <- res_tar[2] / res(r_in)[2]

    # Check if rounding is needed
    if (if_round_fact) {
      fact_x_rounded <- round(fact_x)
      fact_y_rounded <- round(fact_y)
      message(sprintf("  -> Factors: X = %.2f -> %d, Y = %.2f -> %d",
                      fact_x, fact_x_rounded, fact_y, fact_y_rounded))
      fact_x <- fact_x_rounded
      fact_y <- fact_y_rounded
    } else {
      message(sprintf("  -> Factors: X = %.2f, Y = %.2f", fact_x, fact_y))
    }

    message(sprintf("  -> Input raster dimensions: %d x %d",
                    terra::nrow(r_in), terra::ncol(r_in)))

    # Perform aggregation or disaggregation based on factors
    if (fact_x >= 1 && fact_y >= 1) {
      # Downsampling: use aggregate
      r_out <- terra::aggregate(r_in,
                                fact = c(fact_x, fact_y),
                                fun = fun,
                                na.rm = TRUE)
      message(sprintf("  -> Output raster dimensions: %d x %d (after aggregation)",
                      terra::nrow(r_out), terra::ncol(r_out)))
    } else if (fact_x <= 1 && fact_y <= 1) {
      # Upsampling: use disaggregate
      # Calculate disaggregation factors (target resolution -> input resolution)
      disagg_x <- res(r_in)[1] / res_tar[1]
      disagg_y <- res(r_in)[2]  / res_tar[2]

      if (if_round_fact) {
        disagg_x_rounded <- round(disagg_x)
        disagg_y_rounded <- round(disagg_y)
        message(sprintf("  -> Disaggregation factors: X = %.2f -> %d, Y = %.2f -> %d",
                        disagg_x, disagg_x_rounded, disagg_y, disagg_y_rounded))
        disagg_x <- disagg_x_rounded
        disagg_y <- disagg_y_rounded
      } else {
        message(sprintf("  -> Disaggregation factors: X = %.2f, Y = %.2f", disagg_x, disagg_y))
      }

      r_out <- terra::disagg(r_in, fact = c(disagg_x, disagg_y))
      message(sprintf("  -> Output raster dimensions: %d x %d (after disaggregation)",
                      terra::nrow(r_out), terra::ncol(r_out)))
    } else {
      # Mixed case: one factor >1, one <1 - not directly supported o
      warning(sprintf("Mixed factors: X = %.2f, Y = %.2f. ",
                      fact_x, fact_y),
              "Cannot aggregate in one dimension and disaggregate in the other simultaneously. ",
              "Skipping aggregation.")
      r_out <- r_in
    }

  } else {
    # Method 3: No aggregation
    message("No aggregation performed")
    r_out <- r_in
  }

  # ========== Project ==========
  if (if_project && !if_zonal && !is.null(target)) {
    message("Projecting to target")
    r_out <- terra::project(r_out, target, method = "bilinear")
  }

  # ========== Resample ==========
  # Note: resample is only applied when zonal and project is NOT used
  # because zonal already produces exact target grid alignment

  if (if_resample && !if_zonal && !if_project && !is.null(target)) {
    message("Resampling to target grid")
    r_out <- terra::resample(r_out, r_tar, method = "bilinear")
  }

  # ========== Mask ==========
  if (if_mask && !is.null(target)) {
    message("Masking to target extent")
    r_out <- terra::mask(r_out, r_tar)
  }

  # ========== Save ==========
  n_layers <- terra::nlyr(r_out)

  if (!is.null(output)) {
    # If output is a folder path
    if (length(output) == 1 && (dir.exists(output) || !grepl("\\.", basename(output)))) {
      if (!dir.exists(output)) dir.create(output, recursive = TRUE)
      output <- sapply(seq_len(n_layers), function(i) {
        fname <- if (!is.null(names(r_out)[i]) && names(r_out)[i] != "") {
          names(r_out)[i]
        } else {
          paste0("lyr", i)
        }
        file.path(output, paste0(fname, ".nc"))
      })
    }

    # Validate output length
    if (length(output) != n_layers) {
      stop("Length of 'output' must match number of layers in raster.")
    }

    # Save each layer
    for (i in seq_len(n_layers)) {
      lyr <- r_out[[i]]
      out_i <- output[i]
      ext <- tools::file_ext(out_i)

      # Determine variable name
      varname_i <- if (!is.null(names(r_out)[i]) && names(r_out)[i] != "" && names(r_out)[i] != "lyr.1") {
        names(r_out)[i]
      } else {
        if (n_layers == 1) varname else paste0(varname, "_", i)
      }

      # Write based on file extension
      save_raster(lyr, out_i,varname_i)
    }
  }

  # Clean up
  rm(r_in)
  if (exists("r_tar")) rm(r_tar)
  gc()

  if (if_return_raster) return(r_out) else return(output)
}
