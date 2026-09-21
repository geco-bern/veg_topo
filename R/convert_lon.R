convert_lon <- function(x) {
  if (is.numeric(x)) {
    return(paste0(ifelse(x < 0, "W", "E"), sprintf("%03d", abs(x))))
  } else {
    return(ifelse(substr(x, 1, 1) == "W", -1, 1) * as.numeric(substring(x, 2)))
  }
}
