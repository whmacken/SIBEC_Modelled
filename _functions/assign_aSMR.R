#----------------------------------------------------------
# 4. Assign aSMR
#----------------------------------------------------------
# rules_table = fread("./inputsGit/aSMR_Rules_HalfStep_v11_09Dec2019.csv")
# data = sibec.anal

assign_aSMR <- function(
    data,
    rules_table) {
  
  required_rules <- c(
    "SMRLevel",
    "CMD_range",
    "aSMR"
  )
  
  if (!all(required_rules %in% names(rules_table))) {
    stop(
      "rules_table is missing: ",
      paste(
        setdiff(required_rules, names(rules_table)),
        collapse = ", "
      )
    )
  }
  
  if (!"CMD.adj" %in% names(data)) {
    stop("CMD.adj is missing from data.")
  }
  
  data %>%
    mutate(
      SMRLevel = case_when(
        rSMR %in% 0:4 ~ 0,
        rSMR %in% 5:7 ~ rSMR,
        rSMR == 8     ~ NA_real_,
        TRUE          ~ NA_real_
      )
    ) %>%
    
    left_join(
      rules_table,
      by = "SMRLevel"
    ) %>%
    
    filter(
      CMD_range <= CMD.adj |
        rSMR == 8
    ) %>%
    
    group_by(
      PlotNumber,
      rSMR
    ) %>%
    
    slice_max(
      coalesce(CMD_range, -Inf),
      with_ties = FALSE
    ) %>%
    
    ungroup() %>%
    
    mutate(
      aSMR = if_else(
        rSMR == 8,
        0,
        aSMR
      ),
      eda2 = paste0(
        aSMR,
        "/",
        rSNRg
      )
    ) %>%
    
    select(
      -SMRLevel,
      -CMD_range
    )
}


