# wrapper function to get annual total
calc_sw_in <- function(lat, slope = 0, aspect = 0, year = 2020) {
  doy_seq <- 1:(julian_day(year + 1, 1, 1) - julian_day(year, 1, 1))
  daily_rad <- calc_sw_in_daily(lat = lat, slope = slope, aspect = aspect,
                                year = year, doy = doy_seq)
  rowSums(daily_rad) / 1e6  # MJ/m²
  # rowSums(daily_rad)  # J/m²
}


calc_sw_in_daily <- function(
    lat,
    slope = 0.0,
    aspect = 0,
    year = 2001,
    doy
){

  # correction based on SPLASH 2.0
  aspect <- aspect - 180

  ###########################################################################
  # Define constants inside functions to avoid exporting one by one to the cluster
  ###########################################################################
  kA <- 107           # constant for Rl (Monteith & Unsworth, 1990)
  kalb_sw <- 0.17     # shortwave albedo (Federer, 1968)
  kalb_vis <- 0.03    # visible light albedo (Sellers, 1985)
  kb <- 0.20          # constant for Rl (Linacre, 1968; Kramer, 1957)
  kc <- 0.25          # constant for Rs (Linacre, 1968)
  kd <- 0.50          # constant for Rs (Linacre, 1968)
  kfFEC <- 2.04       # from-flux-to-energy, umol/J (Meek et al., 1984)
  kG <- 9.80665       # gravitational acceleration, m/s^2 (Allen, 1973)
  kGsc <- 1360.8      # solar constant, W/m^2 (Kopp & Lean, 2011)
  kL <- 0.0065        # adiabatic lapse rate, K/m (Cavcar, 2000)
  kMa <- 0.028963     # molecular weight of dry air, kg/mol (Tsilingiris, 2008)
  kMv <- 0.01802      # mol. weight of water vapor, kg/mol (Tsilingiris, 2008)
  kSecInDay <- 86400  # number of seconds in a day
  kPo <- 101325       # standard atmosphere, Pa (Allen, 1973)
  kR <- 8.31447       # universal gas constant, J/mol/K (Moldover et al., 1988)
  kTo <- 288.15       # base temperature, K (Berberan-Santos et al., 1997)
  pir <- pi/180       # pi in radians

  # ~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # solar <- list()

  # # obtain orbital parameters
  # orb_out <- orbpar(year)
  #
  # # obliquity
  # keps <- orb_out$obliq
  #
  # # eccentricity
  # ke <- orb_out$eccen
  #
  # # longitude of perihelion
  # komega <- orb_out$long_perihel

  # Paleoclimate variables:
  ke <- 0.01670       # eccentricity of earth's orbit, 2000 CE (Berger 1978)
  keps <- 23.44       # obliquity of earth's elliptic, 2000 CE (Berger 1978)
  komega <- 283       # lon. of perihelion, degrees, 2000 CE (Berger, 1978)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate the number of days in yeark (kN), days
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  kN <- ifelse(year == 0, 365, (julian_day(year + 1, 1, 1) - julian_day(year, 1, 1)))
  # solar$kN <- kN

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 02. Calculate heliocentric longitudes (nu and lambda), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  my_helio <- berger_tls(doy, kN, ke, komega)
  nu <- my_helio$nu
  lam <- my_helio$tls
  # solar$nu_deg <- nu
  # solar$lambda_deg <- lam

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate distance factor (dr), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Berger et al. (1993)
  kee <- ke^2
  rho <- (1 - kee)/(1 + ke*dcos(nu))
  dr <- (1/rho)^2
  # solar$rho <- rho
  # solar$dr <- dr

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate the declination angle (delta), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Woolf (1968)
  delta <- asin(dsin(lam)*dsin(keps))
  delta <- delta/pir
  # solar$delta_deg <- delta

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate variable substitutes (u and v), unitless
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # (Ting) Vectorized computation for multiple latitudes and days

  # Create proper matrices for vectorized operations
  n_days <- length(doy)
  n_lats <- length(lat)

  # Expand all inputs to matrices of correct dimensions
  delta_mat <- matrix(delta, nrow = n_days, ncol = n_lats)
  lat_mat <- matrix(lat, nrow = n_days, ncol = n_lats, byrow = TRUE)
  slope_mat <- matrix(slope, nrow = n_days, ncol = n_lats, byrow = TRUE)
  aspect_mat <- matrix(aspect, nrow = n_days, ncol = n_lats, byrow = TRUE)
  dr_mat <- matrix(dr, nrow = n_days, ncol = n_lats)

  # modification by local slope and aspect
  a <- dsin(delta_mat) * dcos(lat_mat) * dsin(slope_mat) * dcos(aspect_mat) - dsin(delta_mat) * dsin(lat_mat) * dcos(slope_mat)
  b <- dcos(delta_mat) * dcos(lat_mat) * dcos(slope_mat) + dcos(delta_mat) * dsin(lat_mat) * dsin(slope_mat) * dcos(aspect_mat)
  c_val <- dcos(delta_mat) * dsin(slope_mat) * dsin(aspect_mat)  # (Ting) rename

  d <- b^2 + c_val^2 - a^2  # (Ting) use new name
  d[d <= 0] <- 0.000001

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate sin sunset hour angle after Allen, 2006 doi:10.1016/j.agrformet.2006.05.012 0deg is south!!!
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  denominator <- b^2 + c_val^2 # (Ting) use new name
  denominator[denominator == 0] <- 1e-6  # (Ting) Avoid division by zero

  sin_hs <- (a * c_val + b * sqrt(d)) / denominator # use new name
  sin_hs <- pmin(pmax(sin_hs, -1), 1)  # (Ting) Ensure values stay within [-1, 1] to prevent floating-point errors causing NaNs in acos()

  sin_hs[sin_hs < (-1)] <- -1
  sin_hs[sin_hs > (1)] <- 1

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate variable substitutes ru and rv using cos^2(hs)+sin^2(hs) = 1
  # (According to Email David Sandoval 26.05.2025)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ru = -b * sqrt(1.0 - sin_hs^2)
  rv = b

  # Corresponding variable substitute for a flat terrain (slope = 0)
  ru_f <- dsin(delta_mat) * dsin(lat_mat) # (Ting) Vectorized computation for multiple latitudes and days!!
  rv_f <- dcos(delta_mat) * dcos(lat_mat) # (Ting) Vectorized computation for multiple latitudes and days!!


  # (Ting) Vectorized computation for multiple latitudes and days!!
  # correct for anomalous ru, Transparent mountains!
  ru <- ifelse((ru < ru_f) | (ru == 0), ru_f, ru)  # possible issue here????????????????????????????????

  # solar$ru <- ru
  # solar$rv <- rv

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate the sunset hour angle (hs), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note: u/v equals tan(delta) * tan(lat)
  ruv <- ru/rv
  ruv <- pmin(pmax(ruv, -1), 1)  # (Ting) Ensure values stay within [-1, 1] to prevent floating-point errors causing NaNs in acos()

  hs <- acos(-1.0 * ruv) / pir
  hs[ruv >= 1.0] <- 180    # Polar day (no sunset)
  hs[ruv <= -1.0] <- 0     # Polar night (no sunrise)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate daily extraterrestrial radiation (ra_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # ref: Eq. 1.10.3, Duffy & Beckman (1993)
  r_toa <- (kSecInDay/pi) * kGsc * dr_mat * (ru * pir * hs + rv * dsin(hs))

  # solar$r_toa <- r_toa

  # according to Email David Sandoval 26.05.2025
  r_toa[r_toa < 0] <- 0

  # (Ting)
  # r_toa is per unit slope-surface area (J m-2 day-1)
  # project slope-surface irradiance to horizontal-equivalent
  r_toa_horiz_proj <- r_toa / dcos(slope_mat)

  # (Ting) Vectorized computation for multiple latitudes and days!!
  return(t(r_toa_horiz_proj)) # (Ting) change r_toa to r_toa_horiz_proj

}
