# ------calculate_correlations_bywin------------------------------------

#' Calculate windowed correlation statistics
#'
#' Computes correlation statistics for nested spatial windows.
#'
#' @param window_data Data frame from create_spatial_windows()
#' @param x Character, name of first variable (default "twi")
#' @param y Character, name of second variable (default "vegh")
#' @param if_nobs Logical, include observation counts? (default TRUE)
#' @param if_pval Logical, include p-values? (default TRUE)
#' @param if_data Logical, keep raw data? (default FALSE)
#' @param if_peak Logical, calculate peak relationships? (default FALSE)
#' @return A data frame with correlation statistics by window
calculate_correlation_bywin <- function(df_win,
                                        x = "twi",
                                        y = "vegh",
                                        if_nobs = TRUE,
                                        if_pval = TRUE,
                                        if_data = FALSE,
                                        if_peak = FALSE) {
  # Generate full window midpoints (from df_win)
  full_grid <- df_win %>%
    dplyr::select(lon_mid, lat_mid) %>%
    distinct()

  # Compute correlations for actual windows
  df_cor <- df_win |>

    group_by(lon_mid, lat_mid) |>
    tidyr::nest() |>

    mutate(
      # Perform statistical computations
      stats = purrr::map(data, ~{
        df <- .x
        n_obs <- sum(!is.na(df[[x]]) & !is.na(df[[y]]))  # Count only valid observations

        # Initialize result list
        result <- list(
          correlation = NA_real_,
          n_obs = if(if_nobs) n_obs else NULL,
          cor_pval = if(if_pval) NA_real_ else NULL,
          peak = if(if_peak) NA_real_ else NULL
        )

        # Only calculate correlation if there are enough valid observations with variation
        if (n_obs >= 30 && sd(df[[x]], na.rm = TRUE) > 0 && sd(df[[y]], na.rm = TRUE) > 0) {
          test <- cor.test(df[[x]], df[[y]], use = "complete.obs")  # Pearson correlation test; ensure NA-safe
          result$correlation <- test$estimate  # Extract correlation coefficient
          if(if_pval) result$cor_pval <- test$p.value  # Extract p-value
          if(if_peak) result$peak <- identify_peak(df) # Check for peak
        }

        return(result[!sapply(result, is.null)])  # Remove NULL elements
      }),

      # Extract individual fields from the stats list-column
      correlation = purrr::map_dbl(stats, "correlation"),
      n_obs = if(if_nobs) purrr::map_int(stats, "n_obs") else NULL,
      cor_pval = if(if_pval) purrr::map_dbl(stats, "cor_pval") else NULL,
      peak = if(if_peak) purrr::map_dbl(stats, "peak") else NULL
    ) |>
    # Remove NULL columns (those not requested)
    select(-stats) |>
    # Remove columns based on function arguments
    {function(.) {
      if (!if_data) . <- select(., -data)
      # if (!if_nobs) . <- select(., -n_obs)
      # if (!if_pval) . <- select(., -cor_pval)
      # if (!if_peak) . <- select(., -peak)
      .
    }}() |>
    ungroup()

  # Left join to the full grid to fill in missing windows
  df_cor_full <- full_grid |>
    dplyr::left_join(df_cor, by = c("lon_mid", "lat_mid"))

  return(df_cor_full)
}


