
values_per_layers_mean_TEMP <- list()
TEMP <- list()

for(i in 1:5){

TEMP_profile_data <- TEMP_data %>%
  filter(cluster == i) %>%
  left_join(zp_data %>% select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
  left_join(mld_data %>% select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
  left_join(plotData_SIG %>% select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
  mutate(depths_below_zp = NA)
TEMP_profile_data$zp <- round(TEMP_profile_data$zp)


values <- TEMP_profile_data$value

depths = TEMP_profile_data$PRES
dates = TEMP_profile_data$jDay
month = TEMP_profile_data$month
zp = TEMP_profile_data$zp
mld = TEMP_profile_data$mld
depth_bin = 100

depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
# depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
# 
#depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 1000, sep = "_"), depth_seq)
depth_layers <- c("whole_column","productive_layer","depth_0_100","depth_100_200","depth_200_300",
                  "depth_300_400","depth_400_500","depth_500_600","depth_600_700","depth_700_800" )
depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")

data_list_TEMP <- data.frame(values, depths, dates, zp, mld) %>% 
  mutate(depths_below_zp = depths-zp) %>% 
  split(f = .$dates)

values_per_layers_TEMP <- data_list_TEMP %>% ldply(function(df) {
  
  integrated_values_TEMP <- laply(depth_layers, function(depth_layer) {
    if (depth_layer == "whole_column") {index <- which(between(df$depths_below_zp, -Inf, Inf))}
    if (depth_layer == "productive_layer") {index <- which(between(df$depths_below_zp, -Inf, 0))}
    if (grepl("depth", depth_layer)) {
      min_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,1]
      max_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,2]
      index <- which(between(df$depths_below_zp, min_depth_to_keep, max_depth_to_keep))
    }
    mean(df$values[index])
  })
  
  data.frame(df[1,], layer_name = depth_layers, DEPTH = depth_seq, I = integrated_values_TEMP, row.names = NULL)
  
}, .parallel = FALSE, .id = NULL) %>% 
  dplyr::group_by(layer_name) %>%
  dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T),
                Ez_value = c( NA, diff(smoothed_I) / {diff(dates)}),
                Ez_values_smth = smth(Ez_value, window = window, alpha = 2.5, method = "gaussian", tails = T),
                month = ((dates-min(dates))*12)/(max(dates)-min(dates))) %>% 
  ungroup()


  values_per_layers_TEMP <- filter(values_per_layers_TEMP,  DEPTH == "0" |DEPTH == "100" |DEPTH == "200" |DEPTH == "300" |DEPTH == "400" |DEPTH == "500"|DEPTH == "600"|DEPTH == "700" )
  TEMP[[i]] <- print(ggplot(data = values_per_layers_TEMP, aes(x=month, y=smoothed_I, color = DEPTH)) +
                     geom_point() +
                     scale_x_continuous(name = "relative month", expand = c(0,0), limits = c(0,12)) +
                     scale_y_continuous(name = expression(paste(integrated~O[2]~variation~"  [",µmol~O[2]~m^{-2},"] ")), expand = c(0.02,0.02))+
                     theme_bw() +
                     scale_color_manual(name = element_blank(), breaks = c("0", "100", "200", "300", "400","500","600","700"), values=colors,
                                        labels = c("0", "100", "200", "300", "400","500","600","700")))


values_per_layers_mean_TEMP[[i]] <- values_per_layers_TEMP %>%
  group_by(DEPTH) %>%
  summarise(moyenne = mean(I),
            sd = sd(I),
            DEPTH = unique(DEPTH))

}

Temp <- bind_rows(values_per_layers_mean_TEMP[[1]][1,],values_per_layers_mean_TEMP[[2]][1,],
                  values_per_layers_mean_TEMP[[3]][1,],values_per_layers_mean_TEMP[[4]][1,],
                  values_per_layers_mean_TEMP[[5]][1,])

plot(Martin_exp,Temp$moyenne)
