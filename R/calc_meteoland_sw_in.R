calc_meteoland_sw_in <- function(lat_deg, slope_deg, aspect_deg, year = 2020) {

  deg2rad <- function(d) d * pi / 180

  # Convert all inputs to radians once (vectorized)
  lat_rad <- deg2rad(lat_deg)
  slope_rad <- deg2rad(slope_deg)
  aspect_rad <- deg2rad(aspect_deg)

  # Pre-calculate daily solar declination for the entire year (once only)
  dates <- seq(as.Date(paste0(year, "-01-01")),
               as.Date(paste0(year, "-12-31")),
               by = "day")

  date_strings <- format(dates, "%Y-%m-%d")

  # Convert all dates to Julian days at once (vectorized)
  J <- meteoland::radiation_dateStringToJulianDays(date_strings)

  # Pre-compute declination for all days
  delta_values <- sapply(J, meteoland::radiation_solarDeclination)

  # Solar constant (in kW·m-2) - varies daily due to orbital eccentricity
  Sc_values <- sapply(J, meteoland::radiation_solarConstant)

  # Use cumulative summation to avoid storing 365 intermediate values
  n_points <- length(lat_rad)
  potentialRad <- numeric(n_points)

  # Cumulative calculation for each point
  for (i in 1:n_points) {
    daily_sum <- 0
    # Accumulate daily radiation for the entire year
    for (j in seq_along(delta_values)) {
      daily_sum <- daily_sum + meteoland::radiation_potentialRadiation(
        solarConstant = Sc_values[j],
        latrad = lat_rad[i],
        slorad = slope_rad[i],
        asprad = aspect_rad[i],
        delta = delta_values[j]
      )
    }
    # annual potential solar radiation in MJ·m⁻²·year⁻¹
    potentialRad[i] <- daily_sum
  }

  # Project slope-surface irradiance to horizontal-equivalent
  # Avoid division by zero when slope approaches 90°
  cos_slope <- cos(slope_rad)
  cos_slope_safe <- ifelse(abs(cos_slope) < 1e-10, 1e-10, cos_slope)
  potentialRad_proj <- potentialRad / cos_slope_safe

  return(potentialRad_proj)
  # return(potentialRad)

}
