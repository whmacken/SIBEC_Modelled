add_attributes_sibec <- function(
    sibec,
    doy = 172,
    tau = 1,
    download_climate = TRUE,
    vars.choose = grep(
      "^(CMD_an|CMI_an|bFFP|Eref_an|PAS_an|DD5_an|DD18_an|PPT|VPD)",
      list_vars(),
      value = TRUE
    ),
    rules_table = NULL,
    labels = NULL
) {
  
  library(dplyr)
  library(data.table)
  
  #---------------------------------------------------
  # 1. Add SIBEC attributes
  #---------------------------------------------------
  
  sibec <- sibec %>%
    mutate(
      rSNR = gsub("[+-]", "", SNR),
      
      rSNRg = case_when(
        rSNR %in% c("A", "B") ~ "P",
        rSNR %in% c("C") ~ "PM",
        rSNR %in% c("D", "E") ~ "R",
        TRUE ~ "other"
      ),
      
      eda1 = paste0(rSMR, rSNR),
      
      rSMR = as.numeric(gsub("[+-]", "", rSMR)),
      
      rSMR_label = paste0("rSMR", rSMR),
      
      Aspect = as.numeric(Aspect)
    )
  
  #---------------------------------------------------
  # 2. Solar radiation function
  #---------------------------------------------------
  
  calc_solar <- function(data,
                         doy = 172,
                         tau = 1) {
    
    deg2rad <- pi / 180
    
    lat <- data$Latitude
    slope_pct <- data$SlopeGradient
    aspect <- data$Aspect
    
    slope_deg <- atan(slope_pct / 100) * 180 / pi
    
    lat_rad <- lat * deg2rad
    slope_rad <- slope_deg * deg2rad
    aspect_rad <- aspect * deg2rad
    
    decl <- 23.45 * sin(2 * pi * (284 + doy) / 365)
    decl <- decl * deg2rad
    
    alpha <- asin(
      sin(lat_rad) * sin(decl) +
        cos(lat_rad) * cos(decl)
    )
    
    cos_i <-
      sin(alpha) * cos(slope_rad) +
      cos(alpha) * sin(slope_rad) *
      cos(pi - aspect_rad)
    
    cos_i <- pmax(cos_i, 0)
    
    radiation <- 1361 * tau * cos_i
    
    folded_aspect <- abs(180 - abs(aspect - 225))
    
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
        HLI = hli,
        RelativeRadiation = cos_i,
        ClearSkyRadiation = radiation
      )
  }
  
  sibec <- calc_solar(
    sibec,
    doy = doy,
    tau = tau
  )
  
  #---------------------------------------------------
  # 3. Climate data
  #---------------------------------------------------
  
  sibec <- sibec %>%
    arrange(PlotNumber) %>%
    mutate(id = row_number())
  
  if(download_climate) {
    
    my_points <- sibec %>%
      select(
        lat = Latitude,
        lon = Longitude,
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
      return_refperiod = TRUE,
      vars = vars.choose
    )
    
    sibec <- clim_dat %>%
      left_join(sibec, by = "id")
    
    #---------------------------------------------------
    # 4. CMD corrections
    #---------------------------------------------------
    
    sibec$Moist <-
      rowSums(
        sibec[, c(
          "PPT_05",
          "PPT_06",
          "PPT_07",
          "PPT_08",
          "PPT_09"
        )] > 90
      )
    
    sibec <- sibec %>%
      mutate(
        CMD_an = ifelse(
          Moist > 4,
          CMD_an * 0.75,
          CMD_an
        )
      )
    
    sibec <- sibec %>%
      mutate(
        PPT.dorm = PPT_at + PPT_wt,
        Def = pmax(350 - PPT.dorm, 0),
        CMD = CMD_an + Def
      )
    
    sibec <- sibec %>%
      mutate(
        CMD.adj = case_when(
          rSMR == 0 ~ CMD + 500,
          rSMR == 1 ~ CMD + 375,
          rSMR == 2 ~ CMD + 250,
          rSMR == 3 ~ CMD + 125,
          rSMR == 4 ~ CMD,
          rSMR == 5 ~ CMD / 2,
          rSMR == 6 ~ CMD / 4,
          rSMR == 7 ~ CMD / 8,
          rSMR == 8 ~ CMD / 32,
          TRUE ~ NA_real_
        ),
        
        CMD.adj = pmax(CMD.adj, 0)
      )
  }
  
  #---------------------------------------------------
  # 5. Optional aSMR assignment
  #---------------------------------------------------
  
  if(!is.null(rules_table) &&
     !is.null(labels)) {
    
    sibec <- sibec %>%
      left_join(
        rules_table,
        by = c("rSMR" = "SMRLevel")
      ) %>%
      filter(CMD.range <= CMD.adj) %>%
      group_by(PlotNumber, rSMR) %>%
      slice_max(CMD.range) %>%
      ungroup() %>%
      left_join(labels,
                by = "aSMR") %>%
      mutate(
        eda2 = paste0(aSMR_label, "/", rSNRg)
      )
  }
  
  sibec
}