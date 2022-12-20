plotData<-read_csv("data/plotData_climato_DOXY.csv") %>%
  dplyr::mutate(value_avec_offset=value+abs(min(value))) %>%
  dplyr::group_by(cluster,PRES) %>%
  dplyr::mutate(smth_values=smth(value_avec_offset, window = 90, alpha = 2.5, method = "gaussian", tails = T)) %>%
  # filter(cluster==1 & PRES == 1) %>%
  dplyr::mutate(O2_rates = c( NA, diff(smth_values))) %>%
  ungroup()

plotData <- plotData %>%
  group_by(cluster,PRES) %>%
  dplyr::mutate(smth_O2_rates = smth(value, window = 30, alpha = 2.5, method = "gaussian", tails = T)) %>%
  ungroup()

plotData$value = plotData$smth_O2_rates
  

write_csv(plotData,"data/plotData_DOXY_climato_O2_rates.csv")


plotData_temp <- plotData %>%
  filter(cluster==1 & PRES == 400)
 

plotData <- plotData %>%
  mutate(value=O2_rates) %>%
  select(-c(O2_rates,smth_values))
