#----------------------------------------------------------
# 3. Download climate and calculate aCMD
#----------------------------------------------------------

add_climate_data <- function(
    data,
    vars.choose = vars.choose) {
  
  clim_input <-
    data %>%
    arrange(PlotNumber) %>%
    mutate(id = row_number())
  
  my_points <-
    clim_input %>%
    select(
      lat  = Latitude,
      lon  = Longitude,
      elev = Elevation,
      id
    ) %>%
    na.omit() %>%
    as.data.table() %>%
    mutate(
      lon = ifelse(lon > 0, -lon, lon)
    )
  
  clim_dat <- downscale(
    xyz = my_points,
    vars = vars.choose,
    return_refperiod = TRUE
  )
  
  dat <- left_join(
    clim_dat,
    clim_input,
    by = "id"
  )

  dat$Moist <- rowSums(
    dat[, c(
      "PPT_05",
      "PPT_06",
      "PPT_07",
      "PPT_08",
      "PPT_09"
    )] > 90
  )

   dat <- 
    dat %>%
    mutate(
      CMD_an = ifelse(
        Moist > 4,
        CMD_an * 0.75,
        CMD_an
      ),

      PPT.dorm = PPT_at + PPT_wt,

      Def = pmax(
        350 - PPT.dorm,
        0
      ),

      CMD = CMD_an + Def,

      CMD.adj = case_when(
        rSMR == 0 ~ CMD + 500,
        rSMR == 1 ~ CMD + 375,
        rSMR == 2 ~ CMD + 250,
        rSMR == 3 ~ CMD + 125,
        rSMR == 4 ~ CMD,
        rSMR == 5 ~ CMD / 2,
        rSMR == 6 ~ CMD / 4,
        rSMR == 7 ~ CMD / 8,
        rSMR == 8 ~ 0,
        TRUE      ~ NA_real_
      ),

      CMD.site = pmax(CMD.adj, 0)
    )
dat
  #return(dat) 
}
