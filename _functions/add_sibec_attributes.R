#----------------------------------------------------------
# 1. Add SIBEC attributes
#----------------------------------------------------------

add_sibec_attributes <- function(data) {
  
  data %>%
    mutate(
      rSNR = gsub("[+-]", "", SNR),
      
      rSNRg = case_when(
        SNR %in% c("A", "B") ~ "P",
        SNR %in% c("C")      ~ "PM",
        SNR %in% c("D", "E") ~ "R",
        TRUE                            ~ "other"
      ),
      
      eda1 = paste0(rSMR, rSNR),
      
      rSMR = as.numeric(gsub("[+-]", "", rSMR)),
      
      rSMR_label = paste0("rSMR", rSMR),
      
      Aspect = as.numeric(Aspect)
    )
}
