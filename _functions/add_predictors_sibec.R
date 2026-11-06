data = sibec.anal
# add_asmr = TRUE
# rules_table = rules_table

add_predictors_sibec <- function(
    data,
    add_climate = TRUE,
    add_asmr = FALSE,
    rules_table = NULL,
    labels = NULL,
    doy = 172,
    tau = 1) {
  
  data <-
    data %>%
    add_sibec_attributes() %>%
    calc_solar(
      doy = doy,
      tau = tau
    )
  
  if(add_climate) {
    data <- add_climate_data(data)
  }
  
  if(add_asmr) {
    data <- assign_aSMR(
      data,
      rules_table
    )
  }
  
  data
}
