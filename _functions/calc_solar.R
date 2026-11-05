#----------------------------------------------------------
# 2. McCune & Keon heat load / solar radiation
#----------------------------------------------------------

calc_solar <- function(data,
                       doy = 172,
                       tau = 1) {
  
  deg2rad <- pi / 180
  
  slope_deg <- atan(data$SlopeGradient / 100) * 180 / pi
  
  lat_rad    <- data$Latitude * deg2rad
  slope_rad  <- slope_deg * deg2rad
  aspect_rad <- data$Aspect * deg2rad
  
  decl <- 23.45 * sin(2 * pi * (284 + doy) / 365)
  decl <- decl * deg2rad
  
  alpha <- asin(
    sin(lat_rad) * sin(decl) +
      cos(lat_rad) * cos(decl)
  )
  
  cos_i <-
    sin(alpha) * cos(slope_rad) +
    cos(alpha) *
    sin(slope_rad) *
    cos(pi - aspect_rad)
  
  cos_i <- pmax(cos_i, 0)
  
  radiation <- 1361 * tau * cos_i
  
  folded_aspect <- abs(180 - abs(data$Aspect - 225))
  
  hli <-
    0.339 +
    0.808 *
    cos(folded_aspect * deg2rad) *
    cos(slope_rad) -
    0.196 *
    sin(folded_aspect * deg2rad) *
    sin(slope_rad)
  
  hli <- pmin(pmax(hli, 0), 1)
  
  data %>%
    mutate(
      SlopeDeg = slope_deg,
      heat = hli,
      RelativeRadiation = cos_i,
      solar = radiation
    )
}
