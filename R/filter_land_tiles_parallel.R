#' Filter Tiles Containing Land Using a Raster
#'
#' This function takes a data frame of tile extents and a land raster (or raster file path),
#' and returns only the tiles that contain land. The processing is done in parallel, with
#' each worker reading the raster independently to handle large files efficiently.
#'
#' @param tile_grid A data frame with columns xmin, xmax, ymin, ymax defining tile extents.
#' @param target A terra raster object or file path to the raster.
#' @param n_cores Number of parallel cores to use.
#'
#' @return A data frame of tiles from tile_grid that contain land.
#' @export
#'
filter_land_tiles_parallel <- function(tile_grid, target, n_cores = NULL) {
  tile_list <- split(tile_grid, seq_len(nrow(tile_grid)))

  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))  # Ensure cluster is always stopped

  clusterEvalQ(cl, library(terra))

  land_tiles_list <- parLapply(cl, tile_list, function(tile_row, target) {

    land_raster <- if (is.character(target)) rast(target) else target

    tile_extent <- ext(tile_row$xmin, tile_row$xmax, tile_row$ymin, tile_row$ymax)

    cropped_tile <- tryCatch(crop(land_raster, tile_extent), error = function(e) NULL)
    if (is.null(cropped_tile)) return(NULL)

    pixel_values <- values(cropped_tile)
    if (all(is.na(pixel_values))) return(NULL)

    return(tile_row)
  }, target = target)

  land_tiles <- do.call(rbind, Filter(Negate(is.null), land_tiles_list))
  return(land_tiles)
}

