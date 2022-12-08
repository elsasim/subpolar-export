library(plyr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(smoother)
library(tidyr)
library(scales)
library(grid)
library(numDeriv)
library(stats)
library(broom)
library(purrr)

incertitudes_C <- read_csv("data/all/incertitudes_C.csv")
incertitudes_O <- read_csv("data/all/incertitudes_O2.csv")

sec_axis_adjustement_factors <- function(var_to_scale, var_ref, sd_to_scale = NA, sd_ref = NA, log_scale = FALSE) {
  
  if (log_scale) {var_to_scale <- log(var_to_scale)}
  
  index_to_keep <- which(is.finite(var_ref))
  var_ref <- var_ref[index_to_keep]
  sd_ref <- sd_ref[index_to_keep]
  
  index_to_keep <- which(is.finite(var_to_scale))
  var_to_scale <- var_to_scale[index_to_keep]
  sd_to_scale <- sd_to_scale[index_to_keep]
  
  if (any(!is.na(sd_to_scale)) & any(!is.na(sd_ref))) {
    sd_to_scale[which(is.na(sd_to_scale))] <- 0 
    max_var_to_scale <- max(var_to_scale + sd_to_scale, na.rm = T) 
    min_var_to_scale <- min(var_to_scale - sd_to_scale, na.rm = T) 
    sd_ref[which(is.na(sd_ref))] <- 0
    max_var_ref <- max(var_ref + sd_ref, na.rm = T) 
    min_var_ref <- min(var_ref - sd_ref, na.rm = T) 
  } else {
    max_var_to_scale <- max(var_to_scale, na.rm = T) 
    min_var_to_scale <- min(var_to_scale, na.rm = T) 
    max_var_ref <- max(var_ref, na.rm = T) 
    min_var_ref <- min(var_ref, na.rm = T) 
  }
  
  # if (log_scale) {
  # 
  #   # a <- {exp(max_var_ref) - exp(min_var_ref)} / (max_var_to_scale - min_var_to_scale)
  #   # b <- exp(max_var_ref) - a * max_var_to_scale
  #   #
  #   # return(data.frame(diff = a, adjust = b, operation = "scaled var = log(var_to_scale * diff + adjust)",
  #   #                   reverse_operation = "var_to_scale = {exp(scaled_var) - adjust} / diff)"))
  # 
  #   a <- {exp(max_var_ref) - exp(min_var_ref)} / (max_var_to_scale - min_var_to_scale)
  #   b <- exp(max_var_ref) - a * max_var_to_scale
  # 
  #   return(data.frame(diff = a, adjust = b, operation = "scaled var = log(var_to_scale * diff + adjust)",
  #                     reverse_operation = "var_to_scale = {exp(scaled_var) - adjust} / diff)"))
  # }
  
  diff_to_scale <- max_var_to_scale - min_var_to_scale
  diff_to_scale <- ifelse(diff_to_scale == 0, 1 , diff_to_scale)
  diff_ref <- max_var_ref - min_var_ref
  diff <- diff_ref / diff_to_scale
  
  adjust <- (max_var_ref - max_var_to_scale*diff) 
  
  if (log_scale) {
    return(data.frame(diff = diff, adjust = adjust, operation = "scaled var = (log(var_to_scale) * diff) + adjust",
                      trans_axis_operation = "var_to_scale = exp({scaled_var - adjust} / diff)"))
  }
  
  return(data.frame(diff = diff, adjust = adjust, operation = "scaled var = (var_to_scale * diff) + adjust",
                    trans_axis_operation = "var_to_scale = {scaled_var - adjust} / diff)"))
  
}


##############################################################################
####################################################

# Cette fonction permet de determiner les metrics du bloom (apex, onset, climax) a partir des données masquées

mask <- function (profile_data = read_csv("data/all/profile_data_all.csv") %>% arrange(cluster),
                  zp_data = read_csv("data/all/zp_data_all.csv") %>% arrange(cluster),
                  mld_data = read_csv("data/all/mld_data_all.csv") %>% arrange(cluster), 
                  mask_data = read_csv("data/all/mask_data_all.csv") %>% arrange(cluster)%>% filter(PARAM == "POC_masked"),
                  POC_data = read_csv("data/all/POC_data_all.csv") %>% arrange(cluster), 
                  window=30,
                  plot=TRUE){

  Rp_plots <- list()
  
for (i in 1:length(unique(profile_data$cluster))){
  
  print(paste("cluster", i))

####### charger les données
  bbp_profile_data <- profile_data %>%
    filter(PARAM == "BBP700", cluster == i) %>%
    left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>% # zone productive
    left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>% # couche de mélange
    left_join(mask_data %>% dplyr::select(jDay, PRES, POC_mask = value, cluster), by = c("jDay", "PRES", "cluster"))%>% # données masquées
    left_join(POC_data %>% dplyr::select(jDay, PRES, POC = value, cluster), by = c("jDay", "PRES", "cluster")) # données de POC
  
  values <- bbp_profile_data$POC*10^-3 # passer de mg m-3 d-1 à g m-3 d-1 de POC
  
  depths = bbp_profile_data$PRES
  dates = bbp_profile_data$jDay
  month = bbp_profile_data$month
  zp = bbp_profile_data$zp
  mld = bbp_profile_data$mld
  POC_mask = bbp_profile_data$POC_mask*10^-3
  depth_bin = 100 # a adapter en fonction de la taille de couches voulue
  
  
####### analyse 
  
  # intégrer les valeurs masquées
  values_iPOC_mask <- vector() %>% as.numeric()
  for (k in 1:length(unique(dates))){
    range <- which(dates==unique(dates)[k])
    values_iPOC_mask[k] <- sum(POC_mask[range])
  }
  
  #lisser avec une fenetre 'window' definie en argument
  smth_values_iPOC_mask <- smth(values_iPOC_mask, window = window, alpha = 2.5, method = "gaussian", tails = T)
  
  # derivee du stock de POC
  Ez_POC_mask <- c( NA, diff(smth_values_iPOC_mask) / {diff(unique(dates))})
  # derivée sur stock de POC
  rp_values <- Ez_POC_mask/smth_values_iPOC_mask
  
  apex_x <- match(max(smth_values_iPOC_mask),smth_values_iPOC_mask) 
  apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == apex_x)])
  onset_x <- match(min(smth_values_iPOC_mask),smth_values_iPOC_mask)
  onset_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == onset_x)])
  
  range_pos <- vector() # range pos est un vecteur qui contiendra toutes valeurs de croisance entourée de valeurs positive (+/- 5 jours autour)
  
  # climax avec condition qu'au moins 5 jours avant et après aient une croissance positive
  if(onset_x<apex_x){
    for (l in 6:(length(rp_values)-5)){
      if(any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
        range_pos[l] <- NA}
      if(!any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
        range_pos[l] <- rp_values[l]}
    }
    
    climax_x <- match(max(range_pos[onset_x:apex_x], na.rm = T) , range_pos)
    climax_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == climax_x)])
  }
  
  if(onset_x>apex_x){
    for (l in 6:(length(rp_values)-5)){
      if(any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
        range_pos[l] <- NA}
      if(!any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
        range_pos[l] <- rp_values[l]}
    }
    
    climax_x <- match(max(range_pos[-((apex_x+1):(onset_x-1))], na.rm = T) , range_pos)
    climax_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == climax_x)])
    
  }
  
  apex_y <- max(smth_values_iPOC_mask)
  climax_y <- smth_values_iPOC_mask[climax_x]
  onset_y <- smth_values_iPOC_mask[onset_x]
  
  metrics <- data.frame(values = c(onset_x, climax_x, apex_x, onset_y, climax_y, apex_y), 
                        metrics = c("onset_x", "climax_x", "apex_x", "onset_y", "climax_y", "apex_y"),
                        cluster = rep(i,6))
  write_csv(metrics, paste("data/metrics_cluster",i,".csv", sep=""))
  
  ###########
  #####
  scaling_factor_2 <- sec_axis_adjustement_factors(c(0,20), c(-0.02,0.03))
  

  if(plot == TRUE) {
      rp_plot <- data.frame(rp = rp_values, month = unique(month), smth_iPOC_mask = smth_values_iPOC_mask) 
      
      Rp_plots[[i]] <- ggplot(data = rp_plot, aes(x = month)) +
        
        scale_x_continuous(name = "relative month", expand = c(0.01,0.01)) +
        
        scale_y_continuous(expression(paste(r," [",d^{-1},"]")), sec.axis = sec_axis(trans = ~ {. - scaling_factor_2$adjust} / scaling_factor_2$diff , name = expression(paste(iPOC,"  [",g~C~m^{-2},"]")))) + 
        
        geom_line(aes(y = rp), size = 0.5) +
        
        geom_line(aes(y = smth_iPOC_mask * scaling_factor_2$diff + scaling_factor_2$adjust), size = 0.8, color = "#339900") + 
        
        coord_cartesian(xlim = c(0,12), ylim = c(-0.02,0.03)) + 
        
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
        
        geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#339900", size = 1.3) +
        
        geom_vline(xintercept = climax_x_month, linetype = "dotted", color = "#FFCC00", size = 1.3) +
        
        geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF0000", size = 1.3) +
        
        # annotation_custom(grobTree(textGrob(expression(r[max]~"="), x=0.70,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
        # 
        # annotation_custom(grobTree(textGrob(format(rp_plot$rp[climax_x], scientific = T, digits = 2), x=0.81,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
        # 
        # annotation_custom(grobTree(textGrob(paste("onset = j", onset_x), x=0.70,  y=0.86, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
        # annotation_custom(grobTree(textGrob(paste("climax = j", climax_x), x=0.70,  y=0.78, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
        # annotation_custom(grobTree(textGrob(paste("apex = j", apex_x), x=0.70,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
        # annotation_custom(grobTree(textGrob(paste("iPOCmax =", format(max(smth_values_iPOC_mask), scientific = F, digits = 3)), x=0.70,  y=0.62, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
        theme_bw() +
        
        theme(text = element_text(size=25, colour = "black"),
              axis.text.x = element_text(size = 10, colour = "black"), 
              axis.text.y = element_text(size = 10, colour = "black"),
              axis.title.x = element_text(size = 10),
              axis.title.y = element_text(size = 10),
              
              axis.text.y.right = element_text(color = "#339900", size = 10), 
              axis.ticks.y.right = element_line(color ="#339900"),
              axis.line.y.right = element_line(color = "#339900"),
              axis.title.y.right = element_text(color = "#339900", size = 10),
              
              plot.title = element_text(size=12, face="bold.italic", hjust = 0.5)) +
        
        ggtitle(names[i])
  }
  
}
  metric1 <- read_csv("data/metrics_cluster1.csv")
  metric2 <- read_csv("data/metrics_cluster2.csv")
  metric3 <- read_csv("data/metrics_cluster3.csv")
  metric4 <- read_csv("data/metrics_cluster4.csv")
  metric5 <- read_csv("data/metrics_cluster5.csv")
  metric_TOT <- bind_rows(metric1, metric2, metric3, metric4, metric5)
  file.remove("data/metrics_cluster1.csv","data/metrics_cluster2.csv","data/metrics_cluster3.csv","data/metrics_cluster4.csv","data/metrics_cluster5.csv")
  write_csv(metric_TOT,"data/metrics.csv")
  return(Rp_plots)
}

####################################################
# Cette fonction permet de calculer 

iPOC <- function (profile_data, zp_data, mld_data, mask_data, POC_data, window=30){
  
  for (i in 1:length(unique(profile_data$cluster))){
    
    print(paste("cluster", i))
    
####### charger les données
    bbp_profile_data <- profile_data %>%
      filter(PARAM == "BBP700", cluster == i) %>%
      left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>% # zone productive
      left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>% # couche de mélange
      left_join(mask_data %>% dplyr::select(jDay, PRES, POC_mask = value, cluster), by = c("jDay", "PRES", "cluster"))%>% # données masquées
      left_join(POC_data %>% dplyr::select(jDay, PRES, POC = value, cluster), by = c("jDay", "PRES", "cluster")) # données de POC
    
    values <- bbp_profile_data$POC*10^-3 # passer de mg m-3 d-1 à g m-3 d-1 de POC
    
    depths = bbp_profile_data$PRES
    dates = bbp_profile_data$jDay
    month = bbp_profile_data$month
    zp = bbp_profile_data$zp
    mld = bbp_profile_data$mld
    POC_mask = bbp_profile_data$POC_mask*10^-3
    depth_bin = 100 # a adapter en fonction de la taille de couches voulue
    
####### analyse
    depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
    # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
    # 
    #depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 1000, sep = "_"), depth_seq)
    depth_layers <- c("whole_column","productive_layer","depth_0_100","depth_100_200","depth_200_300","depth_300_400","depth_400_500","depth_500_600","depth_600_700","depth_700_800" )
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
        (df$values[index] * {c( diff(df$depths[index]), diff(df$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T) # The integration
      })
      
      data.frame(df[1,], layer_name = depth_layers, DEPTH = depth_seq, I = integrated_values, row.names = NULL)
      
    }, .parallel = FALSE, .id = NULL) %>% 
      dplyr::group_by(layer_name) %>%
      dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T),
                    Ez_value = c( NA, diff(smoothed_I) / {diff(dates)}),
                    Ez_values_smth = smth(Ez_value, window = window, alpha = 2.5, method = "gaussian", tails = T),
                    month = ((dates-min(dates))*12)/(max(dates)-min(dates))) %>% 
      ungroup()
    
    depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
    colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")
    
    # trace bloom phenology with iPOC 0_1000
    I_values_per_layers_ready_to_plot <- values_per_layers %>% pivot_wider(id_cols = c(dates,month, zp, mld), 
                                                                           names_from = "layer_name",values_from = "smoothed_I")
    
    I_plots[[i]] <- ggplot(data = I_values_per_layers_ready_to_plot, aes(x = month))
    
    
    # geom_rect(data = flux_dates_df, aes(xmin = min_date, xmax = max_date, ymin = min_y, ymax = max_y), fill = "green3", alpha = 0.25)
    # geom_area(aes(y = whole_column), fill = "yellow", alpha = 0.5, color = "black", size = 0.15)
    
    
    for (j in 1:8 ) {
      I_plots[[i]] <- I_plots[[i]] + geom_area(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), fill = colors[j], color = colors[j]) +
        geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "black", size = 0.15) +
        geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "white", linetype = "dashed", size = 0.15)
    }
    
    scaling_factor_3 <- sec_axis_adjustement_factors(-c(0,600), c(1,9))
    
    I_plots[[i]] <- I_plots[[i]] +
      geom_line(aes(y = productive_layer, color = "y1"), size=0.5) + 
      scale_x_continuous(name = "relative month", expand = c(0.005,0)) +
      scale_y_continuous(expression(paste(iPOC,"  [",g~C~m^{-2},"]")), breaks =  c(1,3,5,7,9), sec.axis = sec_axis(trans = ~ {. - scaling_factor_3$adjust} / scaling_factor_3$diff *-1, name = "MLD [m]")) + 
      geom_line(aes(y = -mld * scaling_factor_3$diff + scaling_factor_3$adjust, color = "y2"), size = 0.5) + 
      scale_color_manual(values = c("y1" = "#CC0000", "y2" = "#3399FF")) +
      coord_cartesian(xlim = c(0,12), ylim = c(1,9)) + 
      #geom_point(aes(x = apex_x, y = apex_y), shape=21, fill = "#FF3300", size=2, color="#FF3300") +
      geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF0000", size = 1.3) +
      #geom_point(aes(x = onset_x, y = onset_y), shape=21, fill = "#FFCC33", size=2, color="#FFCC33") +
      geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#339900", size = 1.3) +
      #geom_point(aes(x = climax_x, y = climax_y), shape=21, fill = "#FF9900", size=2, color="#FF9900") +
      geom_vline(xintercept = climax_x_month, linetype = "dotted", color = "#FFCC00", size = 1.3) +
      # annotation_custom(grobTree(textGrob(expression(seasonal~E[0]~"="), x=0.63,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) + 
      # annotation_custom(grobTree(textGrob(format(slope_0_100, scientific = T, digits = 2), x=0.84,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) + 
      guides(colour = FALSE) +
      
      theme(text = element_text(size=8, colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(linetype = "solid", fill = NA),
            axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 10, colour = "black"),
            axis.text.x=element_text(angle=0),
            
            axis.text.y.left = element_text(color = "black"), 
            axis.ticks.y.left = element_line(color = "black"),
            axis.line.y.left = element_line(color = "black"),
            axis.title.y.left = element_text(color = "black"),
            
            axis.text.y.right = element_text(color = "#3399FF"), 
            axis.ticks.y.right = element_line(color = "#3399FF"),
            axis.line.y.right = element_line(color = "#3399FF"),
            axis.title.y.right = element_text(color = "#3399FF"),
            plot.title = element_text(size=12, face="bold.italic", hjust = 0.5)) +
      ggtitle(names[i])
    
    
    I_plots[[i]] <- I_plots[[i]] +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_0_100_x)], 
                   y = min_0_100_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_0_100_x)], 
                   yend = max_0_100_y, 
                   color = "#66CC00", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_100_200_x)], 
                   y = min_100_200_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_100_200_x)], 
                   yend = max_100_200_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_200_300_x)], 
                   y = min_200_300_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_200_300_x)], 
                   yend = max_200_300_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_300_400_x)], 
                   y = min_300_400_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_300_400_x)], 
                   yend = max_300_400_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_400_500_x)], 
                   y = min_400_500_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_400_500_x)], 
                   yend = max_400_500_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_500_600_x)], 
                   y = min_500_600_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_500_600_x)], 
                   yend = max_500_600_y, 
                   color = "#339900", size = 0.3)
    
    if(i!=5){
      I_plots[[i]] <- I_plots[[i]] +
        geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_600_700_x)], 
                     y = min_600_700_y, 
                     xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_600_700_x)], 
                     yend = max_600_700_y, 
                     color = "#339900", size = 0.3)+
        geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_700_800_x)], 
                     y = min_700_800_y, 
                     xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_700_800_x)], 
                     yend = max_700_800_y, 
                     color = "#339900", size = 0.3)
    }
    
    export <-  max(I_values_per_layers_ready_to_plot$depth_0_100) /max(I_values_per_layers_ready_to_plot$productive_layer)
    print(export)
  }
}

######################################################################################

attenuation <-  function(profile_data = read_csv("data/all/profile_data_all.csv") %>% arrange(cluster),
                         zp_data = read_csv("data/all/zp_data_all.csv") %>% arrange(cluster),
                         mld_data = read_csv("data/all/mld_data_all.csv") %>% arrange(cluster), 
                         mask_data = read_csv("data/all/mask_data_all.csv") %>% arrange(cluster)%>% filter(PARAM == "POC_masked"),
                         POC_data = read_csv("data/all/POC_data_all.csv") %>% arrange(cluster), 
                         window=30, 
                         metrics = read_csv("data/metrics.csv"),
                         incertitudes_C = read_csv("data/all/incertitudes_C.csv")){
  
  Att_plots <- list()
  Martin_exp <- vector()
  depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
  # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
  # 
  #depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 1000, sep = "_"), depth_seq)
  depth_layers <- c("whole_column","productive_layer","depth_0_100","depth_100_200","depth_200_300","depth_300_400","depth_400_500","depth_500_600","depth_600_700","depth_700_800" )
  # deepest_depth_layer <- depth_layers[which(depth_seq == max_depth_to_consider)]
  # shallowest_depth_layer <- depth_layers[which(depth_seq == 0)]
  
  depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
  colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")
  names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
  
  for (i in 1:length(unique(profile_data$cluster))){
    
    print(paste("cluster", i))
    
####### charger les données
    bbp_profile_data <- profile_data %>%
      filter(PARAM == "BBP700", cluster == i) %>%
      left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>% # zone productive
      left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>% # couche de mélange
      left_join(mask_data %>% dplyr::select(jDay, PRES, POC_mask = value, cluster), by = c("jDay", "PRES", "cluster"))%>% # données masquées
      left_join(POC_data %>% dplyr::select(jDay, PRES, POC = value, cluster), by = c("jDay", "PRES", "cluster")) # données de POC
    
    values <- bbp_profile_data$POC*10^-3 # passer de mg m-3 d-1 à g m-3 d-1 de POC
    
    depths = bbp_profile_data$PRES
    dates = bbp_profile_data$jDay
    month = bbp_profile_data$month
    zp = bbp_profile_data$zp
    mld = bbp_profile_data$mld
    POC_mask = bbp_profile_data$POC_mask*10^-3
    depth_bin = 100 # a adapter en fonction de la taille de couches voulue
    
    metric_temp <- metrics %>% filter(cluster==i)
    onset_x <- metrics$values[1]
    onset_y <- metrics$values[4]
    climax_x <- metrics$values[2]
    climax_y <- metrics$values[5]
    apex_x <- metrics$values[2]
    apex_y <- metrics$values[6]
    
########### analyse
    
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
        (df$values[index] * {c( diff(df$depths[index]), diff(df$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T) # The integration
      })
      
      data.frame(df[1,], layer_name = depth_layers, DEPTH = depth_seq, I = integrated_values, row.names = NULL)
      
    }, .parallel = FALSE, .id = NULL) %>% 
      dplyr::group_by(layer_name) %>%
      dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T),
                    Ez_value = c( NA, diff(smoothed_I) / {diff(dates)}),
                    Ez_values_smth = smth(Ez_value, window = window, alpha = 2.5, method = "gaussian", tails = T),
                    month = ((dates-min(dates))*12)/(max(dates)-min(dates))) %>% 
      ungroup()
    
    # trace bloom phenology with iPOC 0_1000
    whole_column <- values_per_layers %>%
      filter(DEPTH == "whole_column")
    
    iPOC_whole_column <- whole_column$smoothed_I
    
    depth_0_100 <- values_per_layers %>%
      filter(DEPTH == "0")
    iPOC_0_100 <- depth_0_100$smoothed_I
    
    depth_100_200 <- values_per_layers %>%
      filter(DEPTH == "100")
    iPOC_100_200 <- depth_100_200$smoothed_I
    
    depth_200_300 <- values_per_layers %>%
      filter(DEPTH == "200")
    iPOC_200_300 <- depth_200_300$smoothed_I
    
    depth_300_400 <- values_per_layers %>%
      filter(DEPTH == "300")
    iPOC_300_400 <- depth_300_400$smoothed_I
    
    depth_400_500 <- values_per_layers %>%
      filter(DEPTH == "400")
    iPOC_400_500 <- depth_400_500$smoothed_I
    
    depth_500_600 <- values_per_layers %>%
      filter(DEPTH == "500")
    iPOC_500_600 <- depth_500_600$smoothed_I
    
    depth_600_700 <- values_per_layers %>%
      filter(DEPTH == "600")
    iPOC_600_700 <- depth_600_700$smoothed_I
    
    depth_700_800 <- values_per_layers %>%
      filter(DEPTH == "700")
    iPOC_700_800 <- depth_700_800$smoothed_I
    
    
    ######### 0
    min_0_100_x <- match(min(iPOC_0_100[0:apex_x]),iPOC_0_100[0:apex_x])
    min_0_100_y <- iPOC_0_100[min_0_100_x]
    max_0_100_x <- match(max(iPOC_0_100),iPOC_0_100)
    max_0_100_y <- iPOC_0_100[max_0_100_x]
    slope_0_100 <- (max_0_100_y - min_0_100_y) / (max_0_100_x - min_0_100_x) 
    ######### 150
    min_100_200_x <- match(min(iPOC_100_200[0:apex_x]),iPOC_100_200[0:apex_x])
    min_100_200_y <- iPOC_100_200[min_100_200_x]
    max_100_200_x <- match(max(iPOC_100_200),iPOC_100_200)
    max_100_200_y <- iPOC_100_200[max_100_200_x]
    slope_100_200 <- (max_100_200_y - min_100_200_y) / (max_100_200_x - min_100_200_x)
    ######### 300
    min_200_300_x <- match(min(iPOC_200_300[0:apex_x]),iPOC_200_300[0:apex_x])
    min_200_300_y <- iPOC_200_300[min_200_300_x]
    max_200_300_x <- match(max(iPOC_200_300),iPOC_200_300)
    max_200_300_y <- iPOC_200_300[max_200_300_x]
    slope_200_300 <- (max_200_300_y - min_200_300_y) / (max_200_300_x - min_200_300_x)
    ######### 450
    min_300_400_x <- match(min(iPOC_300_400[0:apex_x]),iPOC_300_400[0:apex_x])
    min_300_400_y <- iPOC_300_400[min_300_400_x]
    max_300_400_x <- match(max(iPOC_300_400),iPOC_300_400)
    max_300_400_y <- iPOC_300_400[max_300_400_x]
    slope_300_400 <- (max_300_400_y - min_300_400_y) / (max_300_400_x - min_300_400_x)
    ######### 600
    min_400_500_x <- match(min(iPOC_400_500[0:apex_x]),iPOC_400_500[0:apex_x])
    min_400_500_y <- iPOC_400_500[min_400_500_x]
    max_400_500_x <- match(max(iPOC_400_500),iPOC_400_500)
    max_400_500_y <- iPOC_400_500[max_400_500_x]
    slope_400_500 <- (max_400_500_y - min_400_500_y) / (max_400_500_x - min_400_500_x)
    ######### 750
    min_500_600_x <- match(min(iPOC_500_600[0:apex_x]),iPOC_500_600[0:apex_x])
    min_500_600_y <- iPOC_500_600[min_500_600_x]
    max_500_600_x <- match(max(iPOC_500_600),iPOC_500_600)
    max_500_600_y <- iPOC_500_600[max_500_600_x]
    slope_500_600 <- (max_500_600_y - min_500_600_y) / (max_500_600_x - min_500_600_x)
    
    min_600_700_x <- match(min(iPOC_600_700[0:apex_x]),iPOC_600_700[0:apex_x])
    min_600_700_y <- iPOC_600_700[min_600_700_x]
    max_600_700_x <- match(max(iPOC_600_700),iPOC_600_700)
    max_600_700_y <- iPOC_600_700[max_600_700_x]
    slope_600_700 <- (max_600_700_y - min_600_700_y) / (max_600_700_x - min_600_700_x)
    
    min_700_800_x <- match(min(iPOC_700_800[0:apex_x]),iPOC_700_800[0:apex_x])
    min_700_800_y <- iPOC_700_800[min_700_800_x]
    max_700_800_x <- match(max(iPOC_700_800),iPOC_700_800)
    max_700_800_y <- iPOC_700_800[max_700_800_x]
    slope_700_800 <- (max_700_800_y - min_700_800_y) / (max_700_800_x - min_700_800_x)
    
    
    
    slopes <- data.frame(flux = c(slope_0_100*10^3,slope_100_200*10^3,
                                  slope_200_300*10^3,slope_300_400*10^3,
                                  slope_400_500*10^3,slope_500_600*10^3,
                                  slope_600_700*10^3,slope_700_800*10^3))%>%
      mutate(i_depths = depth_nb+1)
    
    # slopes <- data.frame(flux = c(slope_0_100,slope_100_200,
    #                               slope_200_300,slope_300_400,
    #                               slope_400_500,slope_500_600,
    #                               slope_600_700,slope_700_800))%>%
    #   mutate(i_depths = depth_nb+1) 
    
    if(i==5){
      slopes <- slopes[-(7:8),]
      slopes <- slopes %>%
        mutate(prof = c(50,150,250,350,450,550),
               cluster = rep(i,6))
      incert <- filter(incertitudes_C, cluster==i) %>%
        mutate(prof = c(50,150,250,350,450,550))
    }
    
    if(i!=5){
      slopes <- slopes %>%
        mutate(prof = c(50,150,250,350,450,550,650,750),
               cluster = rep(i,8))
      incert <- filter(incertitudes_C, cluster==i) %>%
        mutate(prof = c(50,150,250,350,450,550,650,750))
    }
    
    
    incert$slope_min <- incert$slope_min
    incert$slope_max <- incert$slope_max
    
    slopes <- slopes %>% mutate(slope_min = incert$slope_min,
                                slope_max = incert$slope_max,
                                diff = slope_max - flux)
    
    write_csv(slopes,paste("data/slopes_cluster",i,".csv", sep=""))
    #################
    # 
    # model <- lm(log(slopes$slope)~log(slopes$i_depths))
    # b = model$coefficients[[2]]
    # Fz0 = model$coefficients[[1]]
    # # plot(log(slopes$slope)~log(slopes$i_depths))
    # # abline(model)
    # ###############
    # if(i==1){
    #   a = 36.84
    #   b = 0.3646
    # }
    # 
    # if(i==2){
    #   a = 24.32 
    #   b = 0.3286
    # }
    # 
    # if(i==3){
    #   a = 19.71
    #   b = 0.3919
    # }
    # 
    # if(i==4){
    #   a = 9.661
    #   b = 0.1769
    # }
    # 
    # if(i==5){
    #   a = 42.9
    #   b = 0.3934
    # }
    
    
    regression_att <- nls(formula = flux ~ a * (prof)^(-b), data = slopes, start = list(a=20, b= 0.3))
    # # regression_att <- drm(formula = flux ~ a * (i_depths/20)^(-b), data = slopes)
    a = summary(regression_att)$parameters[1]
    b = summary(regression_att)$parameters[2]
    # 
    # slopes <- mutate(slopes,prediction = predict(regression_att,slopes))
    
    z = c(50:800)
    regression_data <- data.frame(z, flux = a*(z)^(-b))
    
    # f <- function(prof, a, b) a * (prof)^(-b)
    # slopes <- slopes[-6,]
    Martin_exp[i] <- b 
    
    if(i!=5){
      Att_plots[[i]] <- print(ggplot(data = slopes, aes(x = prof, y = flux)) +
                                geom_col(col="black",fill = colors) +
                                geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
                                # geom_point(shape=21, fill = colors, size=2) +
                                # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
                                geom_line(data = regression_data, aes(x = z, y = flux), color ="darkred", size = 1) +
                                scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,0)) +
                                scale_y_continuous(name = expression(paste(Seasonal~accumulation,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,45)) +
                                coord_flip() +
                                # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                                theme_bw() +
                                theme(axis.text = element_text(size = 11, colour = "black"),
                                      axis.title = element_text(size = 11, colour = "black"),
                                      axis.text.x=element_text(angle=0),
                                      plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                                annotation_custom(grobTree(textGrob(paste("b = ",format(b, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                
                                ggtitle(names[i]))
    }
    
    if(i==5){
      Att_plots[[i]] <- print(ggplot(data = slopes, aes(x = prof, y = flux)) +
                                geom_col(col="black",fill = colors[-(7:8)]) +
                                geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
                                # geom_point(shape=21, fill = colors[-(7:8)], size=2) +
                                # xlim(c(0,40)) +
                                # geom_function(fun = f, args = list(a=a, b=b), color = "darkred",size = 1) +
                                geom_line(data = regression_data[(1:551),], aes(x = z, y = flux), color ="darkred", size = 1) +
                                scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,0)) +
                                scale_y_continuous(name = expression(paste(Accumulation~saisonnière,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,45)) +
                                coord_flip() +
                                # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                                theme_bw() +
                                theme(axis.text = element_text(size = 11, colour = "black"),
                                      axis.title = element_text(size = 11, colour = "black"),
                                      axis.text.x=element_text(angle=0),
                                      plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                                annotation_custom(grobTree(textGrob(paste("b = ",format(b, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                # annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                
                                ggtitle(names[i]))
    }
    
  }
  slopes1 <- read_csv("data/slopes_cluster1.csv")
  slopes2 <- read_csv("data/slopes_cluster2.csv")
  slopes3 <- read_csv("data/slopes_cluster3.csv")
  slopes4 <- read_csv("data/slopes_cluster4.csv")
  slopes5 <- read_csv("data/slopes_cluster5.csv")
  slopes_TOT <- bind_rows(slopes1, slopes2, slopes3, slopes4, slopes5)
  file.remove("data/slopes_cluster1.csv","data/slopes_cluster2.csv","data/slopes_cluster3.csv","data/slopes_cluster4.csv","data/slopes_cluster5.csv")
  write_csv(slopes_TOT,"data/slopes.csv")
  return(Att_plots)
}


#################################################
oxygene <- function(plotData_DOXY = read_csv("data/all/DOXY_data_all.csv", show_col_types = F),
                    zp_data = read_csv("data/all/zp_data_all.csv", show_col_types = F) %>% arrange(cluster),
                    mld_data = read_csv("data/all/mld_data_all.csv", show_col_types = F) %>% arrange(cluster),
                    plotData_SIG = read_csv("data/all/SIG_data_all.csv", show_col_types = F),
                    window=30,
                    incertitudes_O = read_csv("data/all/incertitudes_O2.csv", show_col_types = F)){
  
  Rsub_plots <- list()
  depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
  # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
  # 
  #depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 1000, sep = "_"), depth_seq)
  depth_layers <- c("whole_column","productive_layer","depth_0_100","depth_100_200","depth_200_300","depth_300_400","depth_400_500","depth_500_600","depth_600_700","depth_700_800" )
  # deepest_depth_layer <- depth_layers[which(depth_seq == max_depth_to_consider)]
  # shallowest_depth_layer <- depth_layers[which(depth_seq == 0)]
  
  depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
  colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")
  names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
  
      for (i in 1:length(unique(plotData_DOXY$cluster))){
        
            print(paste("cluster", i))
        
        DOXY_profile_data <- plotData_DOXY %>%
          filter(PARAM == "DOXY", cluster == i) %>%
          left_join(zp_data %>% select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
          left_join(mld_data %>% select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
          left_join(plotData_SIG %>% select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
          mutate(depths_below_zp = NA)
        DOXY_profile_data$zp <- round(DOXY_profile_data$zp)
        
        
        
        
        depths = DOXY_profile_data$PRES
        dates = DOXY_profile_data$jDay
        month = DOXY_profile_data$month
        zp = DOXY_profile_data$zp
        mld = DOXY_profile_data$mld
        depth_bin = 100
        SIG <- DOXY_profile_data$SIG +1000
        values <- DOXY_profile_data$value * SIG # multiplier par densité pour convertir µmol kg-1 en µmol m-3
        
        data_list_O2 <- data.frame(values, depths, dates, zp, mld) %>% 
          mutate(depths_below_zp = depths-zp) %>% 
          split(f = .$dates)
        
      
        
        values_per_layers_O2 <- data_list_O2 %>% ldply(function(df) {
          
          integrated_values_O2 <- laply(depth_layers, function(depth_layer) {
            if (depth_layer == "whole_column") {index <- which(between(df$depths_below_zp, -Inf, Inf))}
            if (depth_layer == "productive_layer") {index <- which(between(df$depths_below_zp, -Inf, 0))}
            if (grepl("depth", depth_layer)) {
              min_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,1]
              max_depth_to_keep <- str_split(gsub("depth_", "", depth_layer), "_", simplify = T)[,2]
              index <- which(between(df$depths_below_zp, min_depth_to_keep, max_depth_to_keep))
            }
            (df$values[index] * {c( diff(df$depths[index]), diff(df$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T)
          })
          
          data.frame(df[1,], layer_name = depth_layers, DEPTH = depth_seq, I = integrated_values_O2, row.names = NULL)
          
        }, .parallel = FALSE, .id = NULL) %>% 
          dplyr::group_by(layer_name) %>%
          dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T),
                        Ez_value = c( NA, diff(smoothed_I) / {diff(dates)}),
                        Ez_values_smth = smth(Ez_value, window = window, alpha = 2.5, method = "gaussian", tails = T),
                        month = ((dates-min(dates))*12)/(max(dates)-min(dates))) %>% 
          ungroup()
        
        
        # if(i !=5){
        #   values_per_layers_O2_temp <- filter(values_per_layers_O2,  DEPTH == "0" |DEPTH == "100" |DEPTH == "200" |DEPTH == "300" |DEPTH == "400" |DEPTH == "500"|DEPTH == "600"|DEPTH == "700" )
        #   O2[[i]] <- print(ggplot(data = values_per_layers_O2_temp, aes(x=month, y=smoothed_I, color = DEPTH)) +
        #                      geom_point() +
        #                      scale_x_continuous(name = "relative month", expand = c(0,0), limits = c(0,12)) +
        #                      scale_y_continuous(name = expression(paste(integrated~O[2]~variation~"  [",µmol~O[2]~m^{-2},"] ")), expand = c(0.02,0.02))+
        #                      theme_bw() +
        #                      scale_color_manual(name = element_blank(), breaks = c("0", "100", "200", "300", "400","500","600","700"), values=colors,
        #                                         labels = c("0", "100", "200", "300", "400","500","600","700")))
        # }
        # 
        # if(i==5){
        #   values_per_layers_O2_temp <- filter(values_per_layers_O2, DEPTH == "0" |DEPTH == "100" |DEPTH == "200" |DEPTH == "300" |DEPTH == "400" |DEPTH == "500" )
        #   O2[[i]] <- print(ggplot(data = values_per_layers_O2_temp, aes(x=month, y=smoothed_I, color = DEPTH)) +
        #                      geom_point() +
        #                      scale_x_continuous(name = "relative month", expand = c(0,0), limits = c(0,12)) +
        #                      scale_y_continuous(name = expression(paste(integrated~O[2]~variation~"  [",µmol~O[2]~m^{-2},"] ")), expand = c(0.02,0.02))+
        #                      theme_bw() +
        #                      scale_color_manual(name = element_blank(), breaks = c("0", "100", "200", "300", "400","500"), values = colors[-(5:6)],
        #                                         labels = c("0", "100", "200", "300", "400","500")))
        # }
        
        
        # trace bloom phenology with iO2 0_1000
        whole_column <- values_per_layers_O2 %>%
          filter(DEPTH == "whole_column")
        
        iO2_whole_column <- whole_column$smoothed_I
        
        depth_0_100 <- values_per_layers_O2 %>%
          filter(DEPTH == "0")
        iO2_0_100 <- depth_0_100$smoothed_I
        
        depth_100_200 <- values_per_layers_O2 %>%
          filter(DEPTH == "100")
        iO2_100_200 <- depth_100_200$smoothed_I
        
        depth_200_300 <- values_per_layers_O2 %>%
          filter(DEPTH == "200")
        iO2_200_300 <- depth_200_300$smoothed_I
        
        depth_300_400 <- values_per_layers_O2 %>%
          filter(DEPTH == "300")
        iO2_300_400 <- depth_300_400$smoothed_I
        
        depth_400_500 <- values_per_layers_O2 %>%
          filter(DEPTH == "400")
        iO2_400_500 <- depth_400_500$smoothed_I
        
        depth_500_600 <- values_per_layers_O2 %>%
          filter(DEPTH == "500")
        iO2_500_600 <- depth_500_600$smoothed_I
        
        depth_600_700 <- values_per_layers_O2 %>%
          filter(DEPTH == "600")
        iO2_600_700 <- depth_600_700$smoothed_I
        
        depth_700_800 <- values_per_layers_O2 %>%
          filter(DEPTH == "700")
        iO2_700_800 <- depth_700_800$smoothed_I
        
        
        ######### 0
        min_0_100_x <- match(min(iO2_0_100),iO2_0_100)
        max_0_100_x <- match(max(iO2_0_100),iO2_0_100)
        min_0_100_y <- iO2_0_100[min_0_100_x]
        max_0_100_y <- iO2_0_100[max_0_100_x]
        if(min_0_100_x>max_0_100_x){
          slope_O2_0_100 <- (min_0_100_y - max_0_100_y) / (min_0_100_x - max_0_100_x)
        }
        if(min_0_100_x<max_0_100_x){
          slope_O2_0_100 <- (min_0_100_y - max_0_100_y) / ((365-max_0_100_x) + min_0_100_x)
        }
        ######### 150
        min_100_200_x <- match(min(iO2_100_200),iO2_100_200)
        min_100_200_y <- iO2_100_200[min_100_200_x]
        max_100_200_x <- match(max(iO2_100_200),iO2_100_200)
        max_100_200_y <- iO2_100_200[max_100_200_x]
        if(min_100_200_x>max_100_200_x){
          slope_O2_100_200 <- (min_100_200_y - max_100_200_y) / (min_100_200_x - max_100_200_x)
        }
        if(min_100_200_x<max_100_200_x){
          slope_O2_100_200 <- (min_100_200_y - max_100_200_y) / ((365-max_100_200_x) + min_100_200_x)
        }
        ######### 300
        min_200_300_x <- match(min(iO2_200_300),iO2_200_300)
        min_200_300_y <- iO2_200_300[min_200_300_x]
        max_200_300_x <- match(max(iO2_200_300),iO2_200_300)
        max_200_300_y <- iO2_200_300[max_200_300_x]
        if(min_200_300_x>max_200_300_x){
          slope_O2_200_300 <- (min_200_300_y - max_200_300_y) / (min_200_300_x - max_200_300_x)
        }
        if(min_200_300_x<max_200_300_x){
          slope_O2_200_300 <- (min_200_300_y - max_200_300_y) / ((365-max_200_300_x) + min_200_300_x)
        }
        ######### 450
        min_300_400_x <- match(min(iO2_300_400),iO2_300_400)
        min_300_400_y <- iO2_300_400[min_300_400_x]
        max_300_400_x <- match(max(iO2_300_400),iO2_300_400)
        max_300_400_y <- iO2_300_400[max_300_400_x]
        if(min_300_400_x>max_300_400_x){
          slope_O2_300_400 <- (min_300_400_y - max_300_400_y) / (min_300_400_x - max_300_400_x)
        }
        if(min_300_400_x<max_300_400_x){
          slope_O2_300_400 <- (min_300_400_y - max_300_400_y) / ((365-max_300_400_x) + min_300_400_x)
        }
        ######### 600
        min_400_500_x <- match(min(iO2_400_500),iO2_400_500)
        min_400_500_y <- iO2_400_500[min_400_500_x]
        max_400_500_x <- match(max(iO2_400_500),iO2_400_500)
        max_400_500_y <- iO2_400_500[max_400_500_x]
        if(min_400_500_x>max_400_500_x){
          slope_O2_400_500 <- (min_400_500_y - max_400_500_y) / (min_400_500_x - max_400_500_x)
        }
        if(min_400_500_x<max_400_500_x){
          slope_O2_400_500 <- (min_400_500_y - max_400_500_y) / ((365-max_400_500_x) + min_400_500_x)
        }
        ######### 750
        min_500_600_x <- match(min(iO2_500_600),iO2_500_600)
        min_500_600_y <- iO2_500_600[min_500_600_x]
        max_500_600_x <- match(max(iO2_500_600),iO2_500_600)
        max_500_600_y <- iO2_500_600[max_500_600_x]
        if(min_500_600_x>max_500_600_x){
          slope_O2_500_600 <- (min_500_600_y - max_500_600_y) / (min_500_600_x - max_500_600_x)
        }
        if(min_500_600_x<max_500_600_x){
          slope_O2_500_600 <- (min_500_600_y - max_500_600_y) / ((365-max_500_600_x) + min_500_600_x)
        }
        
        
        if(i!=5){
          min_600_700_x <- match(min(iO2_600_700),iO2_600_700)
          min_600_700_y <- iO2_600_700[min_600_700_x]
          max_600_700_x <- match(max(iO2_600_700),iO2_600_700)
          max_600_700_y <- iO2_600_700[max_600_700_x]
          if(min_600_700_x>max_600_700_x){
            slope_O2_600_700 <- (min_600_700_y - max_600_700_y) / (min_600_700_x - max_600_700_x)
          }
          if(min_600_700_x<max_600_700_x){
            slope_O2_600_700 <- (min_600_700_y - max_600_700_y) / ((365-max_600_700_x) + min_600_700_x)
          }
          
          min_700_800_x <- match(min(iO2_700_800),iO2_700_800)
          min_700_800_y <- iO2_700_800[min_700_800_x]
          max_700_800_x <- match(max(iO2_700_800),iO2_700_800)
          max_700_800_y <- iO2_700_800[max_700_800_x]
          if(min_700_800_x>max_700_800_x){
            slope_O2_700_800 <- (min_700_800_y - max_700_800_y) / (min_700_800_x - max_700_800_x)
          }
          if(min_700_800_x<max_700_800_x){
            slope_O2_700_800 <- (min_700_800_y - max_700_800_y) / ((365-max_700_800_x) + min_700_800_x)
          }
          
          slopes_O2 <- data.frame(slope = c(slope_O2_0_100,slope_O2_100_200,slope_O2_200_300,
                                            slope_O2_300_400,slope_O2_400_500,slope_O2_500_600, 
                                            slope_O2_600_700,slope_O2_700_800))
          slopes_O2 <- slopes_O2 %>%mutate(i_depths = depth_nb)  
        }
        
        if(i==5){
          slopes_O2 <- data.frame(slope = c(slope_O2_0_100,slope_O2_100_200,slope_O2_200_300,
                                            slope_O2_300_400,slope_O2_400_500,slope_O2_500_600))%>%
            mutate(i_depths = depth_nb[-(7:8)])
        }
        
        slopes_O2$slope[1] <- -(1/1.45)*slopes_O2$slope[1] *10^-3 *12  # x le quotient respiratoire pour passer de O2 en C, 
        slopes_O2$slope[2] <- -(1/1.45)*slopes_O2$slope[2] *10^-3 *12   # x0.001 pour passer de µmol en mmol,
        slopes_O2$slope[3] <- -(1/1.45)*slopes_O2$slope[3] *10^-3 *12   # x12 pour passer mmol de C en mg de C
        slopes_O2$slope[4] <- -(1/1.45)*slopes_O2$slope[4] *10^-3 *12 
        slopes_O2$slope[5] <- -(1/1.45)*slopes_O2$slope[5] *10^-3 *12 
        slopes_O2$slope[6] <- -(1/1.45)*slopes_O2$slope[6] *10^-3 *12 
        
        if(i!=5){
          slopes_O2$slope[7] <- -(1/1.45)*slopes_O2$slope[7] *10^-3 *12 
          slopes_O2$slope[8] <- -(1/1.45)*slopes_O2$slope[8] *10^-3 *12
        }
        
        if(i!=5){
          slopes_O2 <- slopes_O2 %>%
            mutate(prof = c(50,150,250,350,450,550,650,750),
                   cluster = rep(i,8))
          incert_O <- filter(incertitudes_O, cluster==i) %>%
            mutate(prof = c(50,150,250,350,450,550,650,750))
        }
        
        if(i==5){
          slopes_O2 <- slopes_O2 %>%
            mutate(prof = c(50,150,250,350,450,550),
                   cluster = rep(i,6))
          incert_O <- filter(incertitudes_O, cluster==i) %>%
            mutate(prof = c(50,150,250,350,450,550))
        }
        
        incert_O$slope_min <- incert_O$slope_min
        incert_O$slope_max <- incert_O$slope_max
        
        slopes_O2 <- slopes_O2 %>% mutate(slope_min = incert_O$slope_min,
                                          slope_max = incert_O$slope_max,
                                          diff = slope_max - slope)
        write_csv(slopes_O2, paste("data/slopes_O2_cluster",i,".csv",sep=""))
        
        
        if(i!=5){
          Rsub_plots[[i]] <- ggplot(data = slopes_O2, aes(x = prof, y = slope)) +
            geom_col(col="black", fill = colors) +
            geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
            # geom_point(shape=21, fill = colors, size=2) +
            # geom_errorbar(aes(ymin = Rsub-sd_tot, ymax = Rsub+sd_tot), width = 15) +
            scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,-10)) +
            scale_y_continuous(name = expression(paste(Respiration~de~C~intégrée~nette~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0.01,0.01), limits = c(-30,690)) +
            geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
            coord_flip() +
            theme_bw() +
            theme(plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
            ggtitle(names[i])
        }
        
        if(i==5){
          Rsub_plots[[i]] <- ggplot(data = slopes_O2, aes(x = prof, y = slope)) +
            geom_col(col="black", fill = colors[-(7:8)]) +
            geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
            # geom_point(shape=21, fill = colors[-(7:8)], size=1.2) +
            # geom_errorbar(aes(ymin = Rsub-sd_tot, ymax = Rsub+sd_tot), width = 15) +
            scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,-10)) +
            scale_y_continuous(name = expression(paste(Respiration~de~C~intégrée~nette~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0.01,0.01), limits = c(-30,690)) +
            geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
            coord_flip() +
            theme_bw() +
            theme(plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
            ggtitle(names[i])
        }
        
      }
  slopes_O2_1 <- read_csv("data/slopes_O2_cluster1.csv")
  slopes_O2_2 <- read_csv("data/slopes_O2_cluster2.csv")
  slopes_O2_3 <- read_csv("data/slopes_O2_cluster3.csv")
  slopes_O2_4 <- read_csv("data/slopes_O2_cluster4.csv")
  slopes_O2_5 <- read_csv("data/slopes_O2_cluster5.csv")
  slopes_O2_TOT <- bind_rows(slopes_O2_1, slopes_O2_2, slopes_O2_3, slopes_O2_4, slopes_O2_5)
  file.remove("data/slopes_O2_cluster1.csv","data/slopes_O2_cluster2.csv","data/slopes_O2_cluster3.csv","data/slopes_O2_cluster4.csv","data/slopes_O2_cluster5.csv")
  write_csv(slopes_O2_TOT,"data/slopes_O2.csv")

  return(Rsub_plots)
}
  
  ################################# 
  # comparer les slopes d'accumulation de C et de respiration
  
  comparaison_C_O <-  function(slopes = read_csv("data/slopes.csv"),
                               slopes_O2 = read_csv("data/slopes_O2.csv")
                               ){
   Rsub_Att <- bind_cols(att=slopes$flux, slope_min_C = slopes$slope_min, slope_max_C = slopes$slope_max, 
                         Rsub=slopes_O2$slope, slope_min_O = slopes_O2$slope_min, slope_max_O = slopes_O2$slope_max,
                         prof = slopes$prof,
                         cluster = slopes$cluster)
   
   Rsub_Att <- Rsub_Att %>% 
     mutate(depths = prof - 50)
   
   Rsub_Att_plots <- ggplot(data = Rsub_Att, aes(x = att, y = Rsub, group=cluster)) +
     # geom_point(aes(fill=factor(cluster)),shape=21, size = 2) +
     geom_pointrange(aes(ymin = slope_min_O, ymax = slope_max_O, color=factor(cluster)), size = .2) +
     geom_pointrange(aes(xmin = slope_min_C, xmax = slope_max_C, color=factor(cluster)), size = .2) +
     
     scale_x_continuous(name = expression(paste(Accumulation~saisonnière~" [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0.05,0.05)) +
     #scale_x_continuous(name = expression(paste("Transfert efficiency [%]")), expand = c(0.05,0.05)) +
     scale_y_continuous(name = expression(paste(Respiration~nette~de~C~~"[",mg~C~m^{-2}~d^{-1},"]")), expand = c(0.01,0.01)) +
     
     
     # geom_abline(intercept = regression_Rsub_Att[[1]]$coefficients[[1]], slope = regression_Rsub_Att[[1]]$coefficients[[2]], color = "#CC0000", linetype = "dashed") +
     # annotation_custom(grobTree(textGrob(expression("(1)"~R^2~"="), x=0.75,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[1]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # 
     # geom_abline(intercept = regression_Rsub_Att[[2]]$coefficients[[1]], slope = regression_Rsub_Att[[2]]$coefficients[[2]], color = "#00CC00", linetype = "dashed",color = "#FF6666") +
     # annotation_custom(grobTree(textGrob(expression("(2)"~R^2~"="), x=0.75,  y=0.86, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[2]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.86, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # 
     # geom_abline(intercept = regression_Rsub_Att[[3]]$coefficients[[1]], slope = regression_Rsub_Att[[3]]$coefficients[[2]], color = "#FFCC00", linetype = "dashed") +
     # annotation_custom(grobTree(textGrob(expression("(3)"~R^2~"="), x=0.75,  y=0.82, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[3]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.82, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # 
     # geom_abline(intercept = regression_Rsub_Att[[4]]$coefficients[[1]], slope = regression_Rsub_Att[[4]]$coefficients[[2]], color = "#33CCFF", linetype = "dashed") +
     # annotation_custom(grobTree(textGrob(expression("(4)"~R^2~"="), x=0.75,  y=0.78, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[4]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.78, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # 
     # geom_abline(intercept = regression_Rsub_Att[[5]]$coefficients[[1]], slope = regression_Rsub_Att[[5]]$coefficients[[2]], color ="#0066CC", linetype = "dashed") +
     # annotation_custom(grobTree(textGrob(expression("(5)"~R^2~"="), x=0.75,  y=0.74, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     # annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[5]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.74, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     
     theme_bw() +
     # scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
     #                   labels = element_blank(), legend = element_blank()) +
     scale_color_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")) 
   
   
   Rsub_Att$cluster <- as.factor(Rsub_Att$cluster)
   Rsub_depths_plots <- ggplot(data = Rsub_Att, aes(x = depths, y = Rsub, group=cluster, fill=factor(cluster))) +
     geom_point(shape=21, size = 3) +
     scale_x_reverse(name = paste("couches [m]"), expand = c(0.05,0.05)) +
     #scale_x_continuous(name = expression(paste("Transfert efficiency [%]")), expand = c(0.05,0.05)) +
     scale_y_continuous(name = expression(paste(Respiration~nette~de~C~~"[",mg~C~m^{-2}~d^{-1},"]")), expand = c(0.01,0.01)) +
     
     coord_flip() +
     
     theme_bw() +
     scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                       labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) 
   
   
   
   Rsub_Att <- Rsub_Att%>%
     mutate(ratio = att/Rsub * 100) %>%
     mutate(clust_name = c(rep("(1) Austral HCB",8),rep("(2) Austral MCB",8),rep("(3) STZ",8),rep("(4) PN",8),rep("(5) AN",6))) %>%
     filter(att>0) 
   
   
   # resp_plot <- ggplot(data = Rsub_Att, aes(x = clust_name, y = ratio, fill = clust_name)) +
   #   geom_boxplot(aes(fill = clust_name)) +
   #   theme_bw() +
   #   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(fill = factor(depths))) +
   #   scale_fill_manual(values = c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC",colors)) +
   #   theme(axis.title.x = element_blank(),
   #         axis.text.x = element_blank(),
   #         legend.title = element_blank()) +
   #   labs(y="part du POC dans la reminéralisation totale")
   
   resp_plot <- ggplot(data = Rsub_Att, aes(x = as.numeric(prof), y = ratio, color = clust_name)) +
     geom_point(aes(colour = clust_name)) +
     geom_line(aes(colour = clust_name)) +
     scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,0)) +
     scale_y_continuous(name = paste("intensité relative de l'accumulation de POC saisonnière", "\n", 
                                     "par rapport au flux de respiration [%]"), expand = c(0,0),limits = c(0,60)) +
     coord_flip() +
     theme_bw() +
     scale_color_manual(values = c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC")) +
     theme(#axis.title.x = element_blank(),
       #axis.text.x = element_blank(),
       legend.title = element_blank())

}