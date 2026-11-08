make_boundary_data <- function(
    grid_data,
    bgc_filter,
    tree_spp,
    bgc_all,
    snr = c("A", "C", "E"),
    rsmr = c(1, 3, 4, 5, 7),
    n_per_bgc = 5,
    si50 = 0.5,
    purpose = "cold_boundary",
    template = sibec.anal) {
  
  # ----------------------------------------------------------
  # Select and sample boundary locations
  # ----------------------------------------------------------
  
  boundary_dat <- grid_data %>%
    filter({{ bgc_filter }}) %>%
    select(
      PlotNumber = id,
      BGC,
      Latitude = ycoord,
      Longitude = xcoord,
      Elevation = elev
    ) %>%
    mutate(
      Zone = stringr::str_extract(BGC, "[A-Z]+"),
      SlopeGradient = 0,
      Aspect = 0,
      SI50 = si50,
      purpose = purpose
    ) %>%
    slice_sample(
      n = n_per_bgc,
      by = BGC
    )
  
  # ----------------------------------------------------------
  # Keep only BGCs actually represented in sampled data
  # ----------------------------------------------------------
  
  bgc_use <- intersect(
    unique(boundary_dat$BGC),
    bgc_all
  )
  
  # ----------------------------------------------------------
  # Create species x edatopic combinations
  # ----------------------------------------------------------
  
  new_records <- expand.grid(
    SNR = snr,
    rSMR = rsmr,
    TreeSpp = tree_spp,
    BGC = bgc_use,
    stringsAsFactors = FALSE
  )
  
  # ----------------------------------------------------------
  # Expand sampled locations
  # ----------------------------------------------------------
  
  boundary_dat <- boundary_dat %>%
    inner_join(
      new_records,
      by = "BGC"
    ) %>%
    mutate(
      SS_NoSpace = paste0(
        BGC,
        "/",
        rSMR,
        SNR
      ),
      PlotNumber = as.character(PlotNumber),
      rSMR = as.numeric(rSMR)
    )
  
  # ----------------------------------------------------------
  # Match SIBEC structure
  # ----------------------------------------------------------
  
  boundary_dat %>%
    select(
      all_of(names(template))
    )
}