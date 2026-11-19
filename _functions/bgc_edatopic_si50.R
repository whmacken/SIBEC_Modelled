# data  =sibec
# species = "Sx"
# bgc = "SBSmc2"
# min_cell_n = 2
# plot = TRUE

bgc_edatopic_si50 <- function(
    data,
    species,
    bgc = NULL,
    min_cell_n = 3,
    plot = TRUE
) {
  
  #--------------------------------------------
  # 1. Summarize edatopic coverage
  #--------------------------------------------
  
  coverage <- data %>%
    filter(TreeSpp == species, SI50>0.5, !is.na(rSMR)) %>%
    count(
      BGC,
      rSMR,
      SNR,
      name = "n"
    ) %>%
    group_by(BGC) %>% filter(n >= min_cell_n) |> 
    summarise(
      n_cells = n(),
      n_cells_min = sum(n >= min_cell_n),
      n_records = sum(n),
      .groups = "drop"
    ) %>%
    arrange(
      desc(n_cells_min),
      desc(n_cells),
      desc(n_records)
    )
  
  #--------------------------------------------
  # 2. Select best BGC unless supplied
  #--------------------------------------------
  
  if (is.null(bgc)) {
    
    bgc <- coverage %>%
      slice_head(n = 1) %>%
      pull(BGC)
  }
  
  message(
    "Species: ", species,
    " | BGC: ", bgc
  )
  
  #--------------------------------------------
  # 3. Calculate SI50 statistics by cell
  #--------------------------------------------
  
  stats <- data %>%
    filter(
      TreeSpp == species,
      BGC == bgc,
      SI50 > 0.5
    ) %>%
    group_by(
      rSMR,
      SNR
    ) %>%
    summarise(
      n = sum(!is.na(SI50)),
      
      SI50_mean = mean(
        SI50,
        na.rm = TRUE
      ),
      
      SI50_median = median(
        SI50,
        na.rm = TRUE
      ),
      
      SI50_sd = sd(
        SI50,
        na.rm = TRUE
      ),
      
      SI50_min = min(
        SI50,
        na.rm = TRUE
      ),
      
      SI50_max = max(
        SI50,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |> filter(n >= min_cell_n)
  
  #--------------------------------------------
  # 4. Relative SI50
  #--------------------------------------------
  
  max_si <- max(
    stats$SI50_mean,
    na.rm = TRUE
  )
  
  stats <- stats %>%
    mutate(
      Relative_SI50 =
        100 * SI50_mean / max_si,
      
      label = paste0(
        "Mean ", round(SI50_mean, 1),
        "\n",
        round(Relative_SI50), "%",
        "\n",
        round(SI50_min, 1), "–", round(SI50_max, 1),
        "\n",
        "n=", n
      )
    )
  
  #--------------------------------------------
  # 5. Plot
  #--------------------------------------------
  
  if (plot) {
    
    p <- ggplot(
      stats,
      aes(
        x = SNR,
        y = factor(
          rSMR,
          levels = 9:0
        ),
        fill = Relative_SI50
      )
    ) +
      geom_tile(
        colour = "white"
      ) +
      geom_text(
        aes(label = label),
        size = 3
      ) +
      scale_fill_viridis_c(
        limits = c(0, 100),
        name = "Relative\nSI50 (%)"
      ) +
      labs(
        title = paste(
          species,
          "-",
          bgc
        ),
        subtitle = paste0(
          "Mean SI50 relative to maximum cell mean (",
          round(max_si, 1),
          " m)"
        ),
        x = "Nutrient regime",
        y = "Moisture regime"
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank()
      )
    
  } else {
    
    p <- NULL
  }
  
  #--------------------------------------------
  # 6. Return everything
  #--------------------------------------------
  
  list(
    species = species,
    BGC = bgc,
    coverage = coverage,
    statistics = stats,
    plot = p
  )
}