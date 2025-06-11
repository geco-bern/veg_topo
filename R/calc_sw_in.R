# wrapper function to get annual total
calc_sw_in <- function(...){
  sum(
    unlist(
      lapply(1:365, function(doy){
        calc_sw_in_daily(..., doy)
      })
    )
  )
}

calc_sw_in_daily <- function(
  lat,
  slope = 0.0,
  aspect = NA,
  year = 2001,
  doy,
  return_f_toa_terrain = FALSE
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
  # modification by local slope and aspect
  a <- dsin(delta) * dcos(lat) * dsin(slope) * dcos(aspect) - dsin(delta) * dsin(lat) * dcos(slope)
  b <- dcos(delta) * dcos(lat) * dcos(slope) + dcos(delta) * dsin(lat) * dsin(slope) * dcos(aspect)
  c <- dcos(delta) * dsin(slope) * dsin(aspect)
  d <- b^2 + c^2 - a^2
  d[d <= 0] <- 0.000001

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate sin sunset hour angle after Allen, 2006 doi:10.1016/j.agrformet.2006.05.012 0deg is south!!!
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sin_hs <- (a*c + b*sqrt(d))/(b^2 + c^2)
  sin_hs[sin_hs < (-1)] <- -1
  sin_hs[sin_hs > (1)] <- 1

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate variable substitutes ru and rv using cos^2(hs)+sin^2(hs) = 1
  # (According to Email David Sandoval 26.05.2025)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ru = -b * sqrt(1.0 - sin_hs^2)
  rv = b

  # Corresponding variable substitute for a flat terrain (slope = 0)
  ru_f <- dsin(delta) * dsin(lat)
  rv_f <- dcos(delta) * dcos(lat)

  # correct for anomalous ru, Transparent mountains!
  ru <- ifelse(
    ru < ru_f || ru == 0,
    ru_f,
    ru
  )

  # solar$ru <- ru
  # solar$rv <- rv

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate the sunset hour angle (hs), degrees
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Note: u/v equals tan(delta) * tan(lat)
  ruv <- ru/rv
  hs <- ifelse(
    ruv >= 1.0,
    180, # Polar day (no sunset)
    ifelse(
      ruv <= -1.0,
      0, # Polar night (no sunrise)
      acos(-1.0*ruv) / pir
      )
    )

  # for flat earth
  ruv_f <- ru_f/rv_f

  hs_f <- ifelse(
    ruv_f >= 1.0,
    180, # Polar day (no sunset)
    ifelse(
      ruv_f <= -1.0,
      0, # Polar night (no sunrise)
      acos(-1.0*ruv_f) / pir
    )
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Calculate daily extraterrestrial radiation (ra_d), J/m^2
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # ref: Eq. 1.10.3, Duffy & Beckman (1993)
  r_toa <- (86400/pi) * kGsc * dr * (ru * pir * hs + rv * dsin(hs))
  # solar$r_toa <- r_toa

  # according to Email David Sandoval 26.05.2025
  r_toa[r_toa < 0] <- 0

  # for flat earth
  r_toa_f <- (86400/pi) * kGsc * dr * (ru_f * pir * hs_f + rv_f * dsin(hs_f))
  # solar$r_toa <- r_toa

  # according to Email David Sandoval 26.05.2025
  r_toa_f[r_toa_f < 0] <- 0

  # ratio of local incident radiation considering slope and aspect over the
  # hypothetical incident radiation assuming flat earth
  f_toa_terrain <- r_toa / r_toa_f

  if (return_f_toa_terrain){
    return(f_toa_terrain)
  } else {
    return(r_toa)
  }
}
