convert_lat <- function(x) {
  if (is.numeric(x)) {
    return(paste0(ifelse(x < 0, "S", "N"), sprintf("%02d", abs(x))))
  } else {
    return(ifelse(substr(x, 1, 1) == "S", -1, 1) * as.numeric(substring(x, 2)))
  }
}


