I_plots <- list()



POC_sd <- read_csv("data/all/plotData/plotData_sd_POC_Koestner.csv") %>% arrange(cluster)
AOU_sd <- read_csv("data/all/plotData/plotData_sd_AOU.csv") %>% arrange(cluster)
zp_data <- read_csv("data/all/zp_data_all.csv") %>% arrange(cluster)
mld_data <- read_csv("data/all/mld_data_all.csv") %>% arrange(cluster)
mask_data <- read_csv("data/all/mask_data_all.csv") %>% arrange(cluster)%>% filter(PARAM == "Cphyto_masked")

incertitudes_tot_POC <- function(POC_sd, window=30, depth_bin, i, zp, mld) {
  
  # sd_tot_0_800 <- vector()
  # sd_tot_100_800 <- vector()
  # sd_tot_200_800 <- vector()
  # sd_tot_300_800 <- vector()
  # sd_tot_400_800 <- vector()
  # sd_tot_500_800 <- vector()
  # sd_tot_600_800 <- vector()
  # sd_tot_700_800 <- vector()

    # for( i in 1:5){
    
          print(paste("incertitude cluster", i))
          
          sd_profile_data <- POC_sd %>%
            filter(cluster == i)
            # left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
            # left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
            # left_join(mask_data %>% dplyr::select(jDay, PRES, POC_mask = value, cluster), by = c("jDay", "PRES", "cluster"))

          
          # SIG <- DOXY_profile_data$SIG +1000
          # values <- sd_profile_data$value * SIG
          
          values <- sd_profile_data$value * 10^-3
          
          depths = sd_profile_data$PRES
          dates = sd_profile_data$jDay
          # zp = sd_profile_data$zp
          # mld = sd_profile_data$mld
          # POC_mask = sd_profile_data$POC_mask
      
          # if(i==1){
          #   climax_x <- 82
          #   apex_x <- 189
          #   onset_x <- 354
          # }
          # if(i==2){
          #   climax_x <- 54
          #   apex_x <- 178
          #   onset_x <-4
          # }
          # if(i==3){
          #   climax_x <- 56
          #   apex_x <- 117
          #   onset_x <-351
          # }
          # if(i==4){
          #   climax_x <- 72
          #   apex_x <- 89
          #   onset_x <-348
          # }
          # if(i==5){
          #   climax_x <- 76
          #   apex_x <- 149
          #   onset_x <-1
          # }


          
          depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
          # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
          # 
          depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 800, sep = "_"), depth_seq)
          # depth_layers <- c("whole_column","productive_layer","depth_0_100","depth_100_200","depth_200_300","depth_300_400","depth_400_500","depth_500_600","depth_600_700","depth_700_800" )
          # deepest_depth_layer <- depth_layers[which(depth_seq == max_depth_to_consider)]
          # shallowest_depth_layer <- depth_layers[which(depth_seq == 0)]
          
          data_list <- data.frame(values, depths, dates, zp, mld) %>% 
            mutate(depths_below_zp = depths-zp) %>% 
            split(f = .$dates)
          
          values_per_layers <- data_list %>% ldply(function(df) {
            
            integrated_values <- laply(depth_layers, function(depth_layer) {
              if (depth_layer == "whole_column") {index <- which(between(df$depths_below_zp, -Inf, Inf))}
              if (depth_layer == "productive_layer") {index <- which(between(df$depths_below_zp, -Inf, 0))}
              if (grepl("depth", depth_layer)) {
                min_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,1]
                max_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,2]
                index <- which(between(df$depths_below_zp, min_depth_to_keep, max_depth_to_keep))
              }
              sqrt(sum((values[index])^2)) # The integration
            })
            
            data.frame(df[1,], layer_name = depth_layers, layer = depth_seq, I = integrated_values, row.names = NULL)
            
          }, .parallel = FALSE, .id = NULL) %>% 
            dplyr::group_by(layer_name) %>%
            dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>%
                          # Ez_value = c( NA, diff(smoothed_I) / {diff(dates)}),
                          # Ez_values_smth = smth(Ez_value, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>% 
            ungroup()
          
          # depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
          # colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")
          # 
          # I_values_per_layers_ready_to_plot <- values_per_layers %>% pivot_wider(id_cols = c(dates, zp, mld), 
          #                                                                        names_from = "layer_name",values_from = "smoothed_I")
          
          # 
          # sd_tot_0_800 <- I_values_per_layers_ready_to_plot$depth_0_800 
          # sd_tot_100_800 <- I_values_per_layers_ready_to_plot$depth_100_800
          # sd_tot_200_800 <- I_values_per_layers_ready_to_plot$depth_200_800
          # sd_tot_300_800 <- I_values_per_layers_ready_to_plot$depth_300_800
          # sd_tot_400_800 <- I_values_per_layers_ready_to_plot$depth_400_800
          # sd_tot_500_800 <- I_values_per_layers_ready_to_plot$depth_500_800
          # sd_tot_600_800 <- I_values_per_layers_ready_to_plot$depth_600_800
          # sd_tot_700_800 <- I_values_per_layers_ready_to_plot$depth_700_800
          
         return(values_per_layers)
               
    # }
# 
#   return(list(sd_tot_0_800=sd_tot_0_800,
#               sd_tot_100_800=sd_tot_100_800, 
#               sd_tot_200_800=sd_tot_200_800, 
#               sd_tot_300_800=sd_tot_300_800, 
#               sd_tot_400_800=sd_tot_400_800, 
#               sd_tot_500_800=sd_tot_500_800, 
#               sd_tot_600_800=sd_tot_600_800, 
#               sd_tot_700_800=sd_tot_700_800))
#   
}




incertitudes_tot_AOU <- function(AOU_sd, window=30,i, SIG_data, depth_bin, zp, mld) {
  
  # sd_tot_0_800 <- vector()
  # sd_tot_100_800 <- vector()
  # sd_tot_200_800 <- vector()
  # sd_tot_300_800 <- vector()
  # sd_tot_400_800 <- vector()
  # sd_tot_500_800 <- vector()
  # sd_tot_600_800 <- vector()
  # sd_tot_700_800 <- vector()
  
  # for( i in 1:5){
  
  print(paste("incertitude cluster", i))
  
  sd_profile_data <- AOU_sd %>%
    filter(cluster == i)
    # left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
    # left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
    # left_join(mask_data %>% dplyr::select(jDay, PRES, POC_mask = value, cluster), by = c("jDay", "PRES", "cluster"))
  
  
  SIG <- SIG_data %>% filter(cluster == i) 
  SIG <- SIG$value +1000
  values <- sd_profile_data$value * SIG
  
  depths = sd_profile_data$PRES
  dates = sd_profile_data$jDay
  # zp = sd_profile_data$zp
  # mld = sd_profile_data$mld
  # POC_mask = sd_profile_data$POC_mask
  
  # if(i==1){
  #   climax_x <- 82
  #   apex_x <- 189
  #   onset_x <- 354
  # }
  # if(i==2){
  #   climax_x <- 54
  #   apex_x <- 178
  #   onset_x <-4
  # }
  # if(i==3){
  #   climax_x <- 56
  #   apex_x <- 117
  #   onset_x <-351
  # }
  # if(i==4){
  #   climax_x <- 72
  #   apex_x <- 89
  #   onset_x <-348
  # }
  # if(i==5){
  #   climax_x <- 76
  #   apex_x <- 149
  #   onset_x <-1
  # }
  
  
  
  depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
  # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
  # 
  depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 800, sep = "_"), depth_seq)
  # depth_layers <- c("whole_column","productive_layer","depth_0_100","depth_100_200","depth_200_300","depth_300_400","depth_400_500","depth_500_600","depth_600_700","depth_700_800" )
  # deepest_depth_layer <- depth_layers[which(depth_seq == max_depth_to_consider)]
  # shallowest_depth_layer <- depth_layers[which(depth_seq == 0)]
  
  data_list <- data.frame(values, depths, dates, zp, mld) %>% 
    mutate(depths_below_zp = depths-zp) %>% 
    split(f = .$dates)
  
  values_per_layers <- data_list %>% ldply(function(df) {
    
    integrated_values <- laply(depth_layers, function(depth_layer) {
      if (depth_layer == "whole_column") {index <- which(between(df$depths_below_zp, -Inf, Inf))}
      if (depth_layer == "productive_layer") {index <- which(between(df$depths_below_zp, -Inf, 0))}
      if (grepl("depth", depth_layer)) {
        min_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,1]
        max_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,2]
        index <- which(between(df$depths_below_zp, min_depth_to_keep, max_depth_to_keep))
      }
      sqrt(sum((values[index])^2)) # The integration
    })
    
    data.frame(df[1,], layer_name = depth_layers, layer = depth_seq, I = integrated_values, row.names = NULL)
    
  }, .parallel = FALSE, .id = NULL) %>% 
    dplyr::group_by(layer_name) %>%
    dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>%
    # Ez_value = c( NA, diff(smoothed_I) / {diff(dates)}),
    # Ez_values_smth = smth(Ez_value, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>% 
    ungroup()
  
  # depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
  # colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")
  
  # I_values_per_layers_ready_to_plot <- values_per_layers %>% pivot_wider(id_cols = c(dates, zp, mld), 
  #                                                                        names_from = "layer_name",values_from = "smoothed_I")
  # 
  # 
  # sd_tot_0_800 <- I_values_per_layers_ready_to_plot$depth_0_800 
  # sd_tot_100_800 <- I_values_per_layers_ready_to_plot$depth_100_800
  # sd_tot_200_800 <- I_values_per_layers_ready_to_plot$depth_200_800
  # sd_tot_300_800 <- I_values_per_layers_ready_to_plot$depth_300_800
  # sd_tot_400_800 <- I_values_per_layers_ready_to_plot$depth_400_800
  # sd_tot_500_800 <- I_values_per_layers_ready_to_plot$depth_500_800
  # sd_tot_600_800 <- I_values_per_layers_ready_to_plot$depth_600_800
  # sd_tot_700_800 <- I_values_per_layers_ready_to_plot$depth_700_800
  
  
  
  # }
  
  return(values_per_layers)
  
}

