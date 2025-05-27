###########################################################################
# 02. Define functions
###########################################################################
# ************************************************************************
# Name:     julian_day
# Inputs:   - double, year (y)
#           - double, month (m)
#           - double, day of month (i)
# Returns:  double, Julian day
# Features: This function converts a date in the Gregorian calendar
#           to a Julian day number (i.e., a method of consecutative
#           numbering of days---does not have anything to do with
#           the Julian calendar!)
#           * valid for dates after -4712 January 1 (i.e., jde >= 0)
# Ref:      Eq. 7.1 J. Meeus (1991), Chapter 7 "Julian Day", Astronomical
#             Algorithms
# ************************************************************************
julian_day <- function(y, m, i) {
  if (m <= 2) {
    y <- y - 1
    m <- m + 12
  }
  a <- floor(y/100)
  b <- 2 - a + floor(a/4)

  jde <- floor(365.25*(y + 4716)) + floor(30.6001*(m + 1)) + i + b - 1524.5
  return(jde)
}

# ************************************************************************
# Name:     berger_tls
# Inputs:   - double, day of year (n)
#           - double, days in year (N)
# Returns:  numeric list, true anomaly and true longitude
# Features: Returns true anomaly and true longitude for a given day.
# Depends:  - ke ............. eccentricity of earth's orbit, unitless
#           - komega ......... longitude of perihelion
#  Ref:     Berger, A. L. (1978), Long term variations of daily insolation
#             and quaternary climatic changes, J. Atmos. Sci., 35, 2362-2367.
# ************************************************************************
berger_tls <- function(n, N, ke, komega) {

  pir <- pi/180

  # Variable substitutes:
  xee <- ke^2
  xec <- ke^3
  xse <- sqrt(1 - ke^2)

  # Mean longitude for vernal equinox:
  xlam <- (ke/2.0 + xec/8.0)*(1 + xse)*sin(komega*pir) -
    xee/4.0*(0.5 + xse)*sin(2.0*komega*pir) +
    xec/8.0*(1.0/3.0 + xse)*sin(3.0*komega*pir)

  xlam <- 2.0*xlam/pir

  # Mean longitude for day of year:
  dlamm <- xlam + (n - 80.0)*(360.0/N)

  # Mean anomaly:
  anm <- dlamm - komega
  ranm <- anm*pir

  # True anomaly (uncorrected):
  ranv <- ranm + (2.0*ke - xec/4.0)*sin(ranm) +
    5.0/4.0*xee*sin(2.0*ranm) +
    13.0/12.0*xec*sin(3.0*ranm)
  anv <- ranv/pir

  # True longitude:
  my_tls <- anv + komega

  my_tls[my_tls < 0] <- my_tls[my_tls < 0] + 360
  my_tls[my_tls > 360] <- my_tls[my_tls > 360] - 360
  my_nu  <-  my_tls - komega
  my_nu[my_nu < 0] <- my_nu[my_nu < 0]+ 360

  return (list(nu = my_nu, tls = my_tls))
}


# ************************************************************************
# Name:     dcos
# Inputs:   double (d), angle in degrees
# Returns:  double, cosine of angle
# Features: This function calculates the cosine of an angle (d) given
#           in degrees.
# Depends:  pir
# Ref:      This script is based on the Javascript function written by
#           C Johnson, Theoretical Physicist, Univ of Chicago
#           - 'Equation of Time' URL: http://mb-soft.com/public3/equatime.html
#           - Javascript URL: http://mb-soft.com/believe/txx/astro22.js
# ************************************************************************
dcos <- function(d){
  cos(d*pi/180)
}


# ************************************************************************
# Name:     dsin
# Inputs:   double (d), angle in degrees
# Returns:  double, sine of angle
# Features: This function calculates the sine of an angle (d) given
#           in degrees.
# Depends:  pir
# ************************************************************************
dsin <- function(d) {
  sin(d*pi/180)
}
