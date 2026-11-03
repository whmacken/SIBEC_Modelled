#----------------------------------------------------------
# 4. Assign aSMR
#----------------------------------------------------------

assign_aSMR <- function(data,
                        rules_table,
                        labels) {
  
  colnames(rules_table) <-
    c(
      "SMRLevel",
      "CMD.range",
      "aSMR"
    )
  
  data %>%
    left_join(
      rules_table,
      by = c("rSMR" = "SMRLevel")
    ) %>%
    filter(
      CMD.range <= CMD.adj
    ) %>%
    group_by(
      PlotNumber,
      rSMR
    ) %>%
    slice_max(
      CMD.range,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    left_join(
      labels,
      by = "aSMR"
    ) %>%
    mutate(
      eda2 = paste0(
        aSMR_label,
        "/",
        rSNRg
      )
    )
}
