'fonctions finales pour analyse POC au niveau des regions subpolaires'
library(numDeriv)
library(stats)
library(broom)
library(purrr)
library(raster)
library(tidyverse)

library(plyr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(smoother)
library(tidyr)
library(scales)
library(grid)
require(factoextra)
require(tidyverse)
require(lubridate)
require(zoo)
require(RColorBrewer)
require(cowplot)
require(dplyr)
require(tidyr)
require(readr)
require(stringr)
require(data.table)
require(patchwork)
require(graphics)


metrics <- function(mask_data, window=30){

  Rp_plots <- list()
  names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
  
  for (i in 1:length(unique(mask_data$cluster))){
    print(paste("cluster", i))
  
    #variables
    mask_profile <- mask_data %>% filter(cluster==i)
    depths = mask_profile$PRES
    dates = mask_profile$jDay
    month = mask_profile$month
    Cphyto_mask = mask_profile$value
    
    #analyse
    values_iCphyto_mask <- vector() %>% as.numeric()
    for (k in 1:length(unique(dates))){
      range <- which(dates==unique(dates)[k])
      values_iCphyto_mask[k] <- sum(Cphyto_mask[range]*10^-3)
    }
    smth_values_iCphyto_mask <- smth(values_iCphyto_mask, window = window, alpha = 2.5, method = "gaussian", tails = T)
    
    Ez_Cphyto_mask <- c( NA, diff(smth_values_iCphyto_mask) / {diff(unique(dates))})
    rp_values <- Ez_Cphyto_mask/smth_values_iCphyto_mask
    
    #apex
    apex_x <- match(max(smth_values_iCphyto_mask),smth_values_iCphyto_mask) 
    apex_x_month <- unique(mask_profile$month[which(mask_profile$jDay == apex_x)])
    
    #onset
    onset_x <- match(min(smth_values_iCphyto_mask),smth_values_iCphyto_mask)
    onset_x_month <- unique(mask_profile$month[which(mask_profile$jDay == onset_x)])
    
    #climax
    range_pos <- vector()
    
    if(onset_x<apex_x){
      for (l in 6:(length(rp_values)-5)){
        if(any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
          range_pos[l] <- NA}
        if(!any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
          range_pos[l] <- rp_values[l]}
      }
      climax_x <- match(max(range_pos[onset_x:apex_x], na.rm = T) , range_pos)
      climax_x_month <- unique(mask_profile$month[which(mask_profile$jDay == climax_x)])
    }
    
    if(onset_x>apex_x){
      for (l in 6:(length(rp_values)-5)){
        if(any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
          range_pos[l] <- NA}
        if(!any(rp_values[(l-5):(l+5)]<0,na.rm = T)){
          range_pos[l] <- rp_values[l]}
      }
      climax_x <- match(max(range_pos[0:apex_x], na.rm = T) , range_pos)
      climax_x_month <- unique(mask_profile$month[which(mask_profile$jDay == climax_x)])
      
    }
    
    #valeurs des metrics
    apex_y <- max(smth_values_iCphyto_mask) #max de Cphyto
    climax_y <- smth_values_iCphyto_mask[climax_x] #max de croissance
    onset_y <- min(smth_values_iCphyto_mask[0:apex_x]) #debut de croissance
    
    metrics_temp <- data.frame(timings = c(onset_x, climax_x, apex_x), 
                          valeurs = c(onset_y, climax_y, apex_y), 
                          metric = c("onset","climax","apex"), 
                          cluster = c(rep(i,3)))
    
    if(i==1){metrics <- metrics_temp}
    if(i!=1){metrics <- bind_rows(metrics, metrics_temp)}
    
    #plot les taux d'accumulation (rp)
    scaling_factor_2 <- sec_axis_adjustement_factors(c(0,8), c(-0.02,0.03)) #double axe des ordonnées
    rp_plot <- data.frame(rp = rp_values, month = unique(month), smth_iCphyto_mask = smth_values_iCphyto_mask) 
    
    Rp_plots[[i]] <- ggplot(data = rp_plot, aes(x = month)) +
      
      scale_x_continuous(name = "relative month", expand = c(0.01,0.01)) +
      
      scale_y_continuous(expression(paste(r," [",d^{-1},"]")), sec.axis = sec_axis(trans = ~ {. - scaling_factor_2$adjust} / scaling_factor_2$diff , name = expression(paste(iCphyto,"  [",g~C~m^{-2},"]")))) + 
      
      geom_line(aes(y = rp), size = 0.5) +
      
      geom_line(aes(y = smth_iCphyto_mask * scaling_factor_2$diff + scaling_factor_2$adjust), size = 0.8, color = "#339900") + 
      
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
  write_csv(metrics, "data/metrics_cluster.csv")
  
  layout <- c(
    area(1,1),
    area(1,2),
    area(1,3),
    area(2,1),
    area(2,2)
  )
  
  Rp_plots[[1]] <- Rp_plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                                         axis.title.y.right = element_blank(),axis.text.y.right = element_blank())
  Rp_plots[[2]] <- Rp_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank(),
                                         axis.text.x = element_blank(), axis.text.y.left = element_blank(),
                                         axis.title.y.right = element_blank(),axis.text.y.right = element_blank())
  Rp_plots[[3]] <- Rp_plots[[3]] + theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank(),
                                         axis.title.x = element_blank())
  Rp_plots[[4]] <- Rp_plots[[4]] + theme(axis.title.y.right = element_blank(), axis.text.y.right = element_blank(),
                                         axis.title.x = element_blank())
  Rp_plots[[5]] <- Rp_plots[[5]] + theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank())
  
  finalPlot_Rp <- wrap_plots(Rp_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  ggsave(plot = finalPlot_Rp,"figures/Rp_plots.png", width = 11.5, height = 5, dpi = 300)
  
}

metrics(mask_data = read_csv("data/Chla_BBP_data/plotData.csv") %>% arrange(cluster) %>% filter(PARAM == "Cphyto_masked"),
        window = 30)



POC_flux <- function(profile_data, zp_data, mld_data, metrics, window, POC_sd){
  
  I_plots <- list()
  # Att_plots <- list()
  Martin_exp <- vector()
  names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
  
  for (i in 1:length(unique(profile_data$cluster))){
    
    all_data <- profile_data %>%
      filter(PARAM == "POC_Koestner", cluster == i) %>%
      left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster"))
    
    #values
    values <- all_data$value*10^-3 # mg.m-3 en g.m-3
    depths = all_data$PRES
    dates = all_data$jDay
    month = all_data$month
    zp = all_data$zp
    mld = all_data$mld
    depth_bin = 100
    
    metrics_cluster <- metrics %>% filter(cluster == i)
    
    onset_x <- metrics_cluster$timings[1]
    climax_x <- metrics_cluster$timings[2]
    apex_x <- metrics_cluster$timings[3]
    
    onset_x_month <- unique(profile_data$month[which(profile_data$jDay==onset_x)])
    climax_x_month <- unique(profile_data$month[which(profile_data$jDay==climax_x)])
    apex_x_month <- unique(profile_data$month[which(profile_data$jDay==apex_x)])
    
    onset_y <- metrics_cluster$valeurs[1]
    climax_y <- metrics_cluster$valeurs[2]
    apex_y <- metrics_cluster$valeurs[3]
    
    # definir prof max sous zp pour le cluster
    prof_max <- 1000-max(zp)
    
    if(prof_max>850){
      # analyse
      depth_seq <- c("whole_column", "productive_layer", seq(0, 850-depth_bin, by = depth_bin))
      depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 850, sep = "_"), depth_seq)
    
    }else if(prof_max>800){
    # analyse
    depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
    depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 800, sep = "_"), depth_seq)
    
    }else if(prof_max<800){
      # analyse
      depth_seq <- c("whole_column", "productive_layer", seq(0, 650-depth_bin, by = depth_bin))
      depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 650, sep = "_"), depth_seq)}
    
    

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
    
    #calcul des slopes pour chaque couche
    iPOC <- values_per_layers %>% 
      dplyr::group_by(DEPTH) %>%
      dplyr::summarise(min_y = min(smoothed_I[0:182], na.rm=T),
                       max_y = max(smoothed_I, na.rm=T),
                       min_x = match(min(smoothed_I[0:182], na.rm=T),smoothed_I[0:182]),
                       max_x = match(max(smoothed_I, na.rm=T),smoothed_I),
                       Dprod = if(min_x>max_x){Dprop = 365 + max_x - min_x}
                       else if(min_x<max_x){Dprop = max_x - min_x},
                       slope = if(min_x>max_x){slope = (max_y - min_y) / (365 + max_x - min_x) *10^3}
                       else if(min_x<max_x){slope = (max_y - min_y) / (max_x - min_x)*10^3}) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(slope_annual = slope * Dprod) %>%
      dplyr::filter(DEPTH!= "whole_column" & DEPTH!= "productive_layer") %>% 
      dplyr::arrange(DEPTH)
    
    
    sd_tot <-  incertitudes_tot_POC(POC_sd, window=30, i, depth_bin=depth_bin, zp=zp, mld=mld, depth_seq=depth_seq, depth_layers=depth_layers)%>%
      filter(layer!= "whole_column" & layer!= "productive_layer")
      
    iPOC <- iPOC %>%
      dplyr::group_by(DEPTH) %>%
      dplyr::mutate(incert_min = sd_tot %>% filter(dates == min_x & layer == DEPTH ) %>% dplyr::select(smoothed_I) %>% as.numeric(),
                    incert_max = sd_tot %>% filter(dates == max_x & layer == DEPTH ) %>% dplyr::select(smoothed_I)%>% as.numeric(),
                    min_sd_low = min_y - incert_min,
                    min_sd_high = min_y + incert_min,
                    max_sd_low = max_y - incert_max,
                    max_sd_high = max_y + incert_max,
                    slope_max = (max_sd_high - min_sd_low) / (max_x - min_x) * 10^3,
                    slope_min = (max_sd_low - min_sd_high) / (max_x - min_x) * 10^3) %>%
      ungroup() %>%
      dplyr::mutate(cluster = i,
                    Teff = slope_annual / slope_annual[1])
    
    
    
    # agreger les slopes ensemble a chaque iterration pour avoir un fichier final
    # avec tout
    if(i==1){
      slopes_final <- iPOC
    }
    
    if(i!=1){
      slopes_final <- bind_rows(slopes_final, iPOC)
    }
    
    ###### Iplot ###### 
    I_values_per_layers_ready_to_plot <- values_per_layers %>% pivot_wider(id_cols = c(dates,month, zp, mld), 
                                                                           names_from = "layer_name",values_from = "smoothed_I")
    I_plots[[i]] <- ggplot(data = I_values_per_layers_ready_to_plot, aes(x = month)) +
      geom_area(aes(y = whole_column), fill = "wheat") 
      # geom_line(aes(y = whole_column), col = "grey", linetype = "dashed", size = 0.15)
    
    
    for (j in 1:length(colors)) {
      I_plots[[i]] <- I_plots[[i]] + geom_area(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), fill = colors[j], color = colors[j])
        # geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "black", size = 0.15) +
        # geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "white", linetype = "dashed", size = 0.15)
    }
    
    scaling_factor_3 <- sec_axis_adjustement_factors(-c(0,600), c(1,30))
    
    I_plots[[i]] <- I_plots[[i]] +
      geom_line(aes(y = productive_layer, color = "y1"), size=0.5) +
      scale_x_continuous(name = "relative month", expand = c(0.005,0)) +
      scale_y_continuous(expression(paste(iPOC,"  [",g~C~m^{-2},"]")), breaks =  c(1,5,10,15,20,25,30), sec.axis = sec_axis(trans = ~ {. - scaling_factor_3$adjust} / scaling_factor_3$diff *-1, name = "MLD [m]")) + 
      geom_line(aes(y = -mld * scaling_factor_3$diff + scaling_factor_3$adjust, color = "y2"), size = 0.5) + 
      scale_color_manual(values = c("y1" = "#CC0000", "y2" = "#3399FF")) +
      coord_cartesian(xlim = c(0,12), ylim = c(1,30)) + 
      #geom_point(aes(x = apex_x, y = apex_y), shape=21, fill = "#FF3300", size=2, color="#FF3300") +
      # geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF0000", size = 1.3) +
      #geom_point(aes(x = onset_x, y = onset_y), shape=21, fill = "#FFCC33", size=2, color="#FFCC33") +
      # geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#339900", size = 1.3) +
      #geom_point(aes(x = climax_x, y = climax_y), shape=21, fill = "#FF9900", size=2, color="#FF9900") +
      # geom_vline(xintercept = climax_x_month, linetype = "dotted", color = "#FFCC00", size = 1.3) +
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
    
    
    # I_plots[[i]] <- I_plots[[i]] +
    #   geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_0_800_x)], 
    #                y = min_0_800_y, 
    #                xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_0_800_x)], 
    #                yend = max_0_800_y, 
    #                color = "#66CC00", size = 0.3) +
    #   geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_100_800_x)], 
    #                y = min_100_800_y, 
    #                xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_100_800_x)], 
    #                yend = max_100_800_y, 
    #                color = "#339900", size = 0.3) +
    #   geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_200_800_x)], 
    #                y = min_200_800_y, 
    #                xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_200_800_x)], 
    #                yend = max_200_800_y, 
    #                color = "#339900", size = 0.3) +
    #   geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_300_800_x)], 
    #                y = min_300_800_y, 
    #                xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_300_800_x)], 
    #                yend = max_300_800_y, 
    #                color = "#339900", size = 0.3) +
    #   geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_400_800_x)], 
    #                y = min_400_800_y, 
    #                xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_400_800_x)], 
    #                yend = max_400_800_y, 
    #                color = "#339900", size = 0.3) +
    #   geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_500_800_x)], 
    #                y = min_500_800_y, 
    #                xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_500_800_x)], 
    #                yend = max_500_800_y, 
    #                color = "#339900", size = 0.3)
    # 
    # if(i!=5){
    #   I_plots[[i]] <- I_plots[[i]] +
    #     geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_600_800_x)], 
    #                  y = min_600_800_y, 
    #                  xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_600_800_x)], 
    #                  yend = max_600_800_y, 
    #                  color = "#339900", size = 0.3)+
    #     geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_700_800_x)], 
    #                  y = min_700_800_y, 
    #                  xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_700_800_x)], 
    #                  yend = max_700_800_y, 
    #                  color = "#339900", size = 0.3)
    # }
    

    
    ###### Att plots ###### 
    iPOC <- iPOC %>% transform(DEPTH = as.numeric(DEPTH))
    
    regression_att <- nls(formula = slope ~ slope[1] * (DEPTH+1) ^(-b), data = iPOC, start = list(b=1))
    b = summary(regression_att)$parameters[1]
    
    z = c(0:800)
    if(i==1){ regression_data <- data.frame(z = z, slope = iPOC$slope[1] * (z) ^(-b), cluster = i) }
    if(i!=1){ regression_data_temp <- data.frame(z = z, slope = iPOC$slope[1] * (z) ^(-b), cluster = i)
              regression_data <- bind_rows(regression_data, regression_data_temp)}
    
    Martin_exp[i] <- b 
    
  }
  
  slopes_final$DEPTH <- as.numeric(slopes_final$DEPTH)
  
  

  Att_plots <- ggplot(data = slopes_final, aes(x = DEPTH, y = slope, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
                        geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
                        geom_point() +
                        geom_line() +
                        # geom_point(shape=21, fill = colors, size=2) +
                        # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
                        # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
                        scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
                        scale_y_continuous(name = expression(paste(Carbon~flux,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,90)) +
                        coord_flip() +
                        # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                        theme_bw() +
                        theme(axis.text = element_text(size = 11, colour = "black"),
                              axis.title = element_text(size = 11, colour = "black"),
                              axis.text.x=element_text(angle=0),
                              plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                        scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                                          labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
                        scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  Att_plots_annual <- ggplot(data = slopes_final, aes(x = DEPTH, y = slope_annual*10^-3, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
    geom_point() +
    geom_line() +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
    # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
    scale_y_continuous(name = expression(paste(Annual~carbon~flux,"  [",g~C~m^{-2}~y^{-1},"]")), expand = c(0,0),limits = c(0,10)) +
    coord_flip() +
    # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
    theme_bw() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.text.x=element_text(angle=0),
          plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  Teff_plots_annual <- ggplot(data = slopes_final, aes(x = DEPTH, y = Teff, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
    geom_point() +
    geom_line() +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
    # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0.01,0.01), limits = c(850,0)) +
    scale_y_continuous(name = expression(paste(Teff,"  [%]")), expand = c(0.01,0.01),limits = c(0,1)) +
    coord_flip() +
    # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
    theme_bw() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.text.x=element_text(angle=0),
          plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
                  
   
  
  layout <- c(
    area(1,1),
    area(1,2),
    area(1,3),
    area(2,1),
    area(2,2)
  )
  
  I_plots[[1]] <- I_plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                                       axis.title.y.right = element_blank(),axis.text.y.right = element_blank())
  I_plots[[2]] <- I_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank(),
                                       axis.text.x = element_blank(), axis.text.y.left = element_blank(),
                                       axis.title.y.right = element_blank(),axis.text.y.right = element_blank())
  I_plots[[3]] <- I_plots[[3]] + theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank(),
                                       axis.title.x = element_blank())
  I_plots[[4]] <- I_plots[[4]] + theme(axis.title.y.right = element_blank(), axis.text.y.right = element_blank(),
                                       axis.title.x = element_blank())
  I_plots[[5]] <- I_plots[[5]] + theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank())
  
  finalPlot_I <- wrap_plots(I_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  ggsave(plot = finalPlot_I,"figures/I_plots_POC.png", width = 11.5, height = 5, dpi = 300)
  
  
  ggsave(plot = Att_plots,"figures/Att_plots_POC.png", width = 6, height = 6, dpi = 300)
  write_csv(slopes_final, "data/slopes_POC.csv")
  
  ggsave(plot = Att_plots_annual,"figures/Att_plots_annual_POC.png", width = 6, height = 6, dpi = 300)

}

# charger fonction incertitude_C
POC_flux(profile_data = read_csv("data/Chla_BBP_data/plotData.csv", show_col_types = F) %>% arrange(cluster),
         zp_data = read_csv("data/zp_data_all.csv", show_col_types = F) %>% arrange(cluster),
         mld_data = read_csv("data/mld_data_all.csv", show_col_types = F) %>% arrange(cluster),
         metrics = read_csv("data/metrics_cluster.csv", show_col_types = F) %>% arrange(cluster),
         window = 30,
         POC_sd = read_csv("data/Chla_BBP_data/plotData_sd_POC_Koestner.csv", show_col_types = F) %>% arrange(cluster))



# AOU_flux <- function(data_DOXY, zp_data, mld_data, SIG_data, window, AOU_sd){
#   
#   I_plots <- list()
#   names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
#   
#   for (i in 1:length(unique(data_DOXY$cluster))){
#     
#     DOXY_profile_data <- data_DOXY %>%
#       filter(PARAM == "AOU", cluster == i) %>%
#       left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
#       left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
#       left_join(SIG_data %>% dplyr::select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) 
#       # mutate(depths_below_zp = NA)
#     # DOXY_profile_data$zp <- round(DOXY_profile_data$zp)
#     
#     
#     depths = DOXY_profile_data$PRES
#     dates = DOXY_profile_data$jDay
#     month = DOXY_profile_data$month
#     zp = DOXY_profile_data$zp
#     mld = DOXY_profile_data$mld
#     depth_bin = 50
#     SIG <- DOXY_profile_data$SIG +1000
#     values <- DOXY_profile_data$value * SIG
#     depths_below_zp = round(depths-zp)
#     
#     # definir prof max sous zp pour le cluster
#     prof_max <- 1000-max(zp)
#     
#     if(prof_max>850){
#       # analyse
#       depth_seq <- c("whole_column", "productive_layer", seq(0, 850-depth_bin, by = depth_bin))
#       depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 850, sep = "_"), depth_seq)
#       
#       # depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, (as.numeric(depth_seq)+depth_bin), sep = "_"), depth_seq)
#       
#     }else if(prof_max>800){
#       # analyse
#       depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
#       depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 800, sep = "_"), depth_seq)
#       # depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, (as.numeric(depth_seq)+depth_bin), sep = "_"), depth_seq)
#       
#     }else if(prof_max<800){
#       # analyse
#       depth_seq <- c("whole_column", "productive_layer", seq(0, 650-depth_bin, by = depth_bin))
#       depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 650, sep = "_"), depth_seq)}
#       # depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, (as.numeric(depth_seq)+depth_bin), sep = "_"), depth_seq)}
#     
#       
#     data_list_O2 <- data.frame(values, depths, dates, zp, mld, depths_below_zp) %>% 
#       arrange(depths_below_zp) %>%
#       # filter(between(depths_below_zp,0,as.numeric(last(depth_seq))+50)) %>%
#       split(f = .$depths_below_zp)
#     
#     
#     # following Su et al 2022
#     ZP=round(max(zp))+1
#     maxi=as.numeric(last(depth_seq))+50+ZP
#     resp_each_depth <- data_list_O2[ZP:maxi] %>%
#       llply(function(x){
#       metrics <- data.frame(x) %>%
#       dplyr::mutate(smoothed_values = smth(values, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>% # smoothed data
#       dplyr::summarise(values=values,
#                        dates=dates,
#                        min_y = min(smoothed_values, na.rm=T), # minimum on smoothed data
#                        max_y = max(smoothed_values, na.rm=T), # maximum on smoothed data
#                        min_x = match(min(smoothed_values, na.rm=T),smoothed_values) %>% as.numeric(), # minimum on smoothed data
#                        max_x = match(max(smoothed_values, na.rm=T),smoothed_values) %>% as.numeric(), # maximum on smoothed data
#                        Dprod = if(min_x>max_x){Dprod = 365 + max_x - min_x} # productive period between min and max
#                        else if(min_x<max_x){Dprod = max_x - min_x},
#                        min_supp_max = if(min_x>max_x){min_supp_max = TRUE}
#                        else if(min_x<max_x){min_supp_max = FALSE}) %>%
#       dplyr::mutate(dates_recal = ifelse(min_supp_max==T & between(dates,0,max_x),dates+365,dates), # recal data if max is before min pretending continuous TS
#                       max_x_recal = ifelse(min_supp_max==T & between(dates,0,max_x),max_x+365,max_x)) %>%
#       dplyr::filter(dates_recal>min_x | dates_recal<max_x_recal) %>%
#       dplyr::mutate(slopes = summary(lm(values~ dates_recal, data = .))$coefficients[2]* (1/1.45) *10^-3 *12, # linear regression coeff on unsmoothed data and µmol/m3 O2 to mg C/m3
#                     intercept = summary(lm(values~ dates_recal, data = .))$coefficients[1],
#                     slope_annual = slopes * Dprod)
#       })
#      
#       respiration <- resp_each_depth %>% ldply(function(x){
#         values_to_plot <- data.frame(x) %>%
#           dplyr::summarise(R = unique(slopes),
#                            annual_R = unique(slope_annual))
#       })
#       
#       names(respiration)[1] <- "depths"
#       respiration$depths <- as.numeric(respiration$depths)
#       respiration <- respiration %>% dplyr::mutate(cluster=i)
#       
#       
#     ANCP <- data.frame(ANCP=rep(NA,length(depth_layers)), 
#                        depth=rep(NA,length(depth_layers)),
#                        resp=rep(NA,length(depth_layers)),
#                        cluster=rep(i,length(depth_layers)))
#    for(l in 1:length(depth_layers)){
#         if (grepl("depth", depth_layers[l])) {
#           min_depth_to_keep <- str_split(gsub("depth_", "", depth_layers[l]), "_", simplify = T)[,1] %>% as.numeric()
#           max_depth_to_keep <- str_split(gsub("depth_", "", depth_layers[l]), "_", simplify = T)[,2] %>% as.numeric()
#           index <- which(between(respiration$depths, min_depth_to_keep, max_depth_to_keep))
#           ANCP$resp[l] <- (respiration$R[index] * {c( diff(respiration$depths[index]), diff(respiration$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T) #mgC/m2/jour
#           ANCP$depth[l] <- min_depth_to_keep
#           ANCP$ANCP[l] <- (respiration$annual_R[index] * {c( diff(respiration$depths[index]), diff(respiration$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T)*10^-3 #gC/m2/an
#         }
#       }
#       
#     # agreger les slopes ensemble a chaque iterration pour avoir un fichier final
#     # avec tout
#     
#     if(i==1){
#       slopes_final <- ANCP
#       resp_final <- respiration
#     }
#     
#     if(i!=1){
#       slopes_final <- bind_rows(slopes_final, ANCP)
#       resp_final <- bind_rows(resp_final, respiration)
#     }
#     
#   }
#   
#   resp_each_depth_plot <- ggplot(data = resp_final, aes(x = depths, y = R, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
#     geom_line() +
#     # geom_errorbar(aes(ymin = estimate-sd, ymax = estimate+sd), width = 15) +
#     scale_x_reverse(name = "Pressure under Zp [dbar]", expand = c(0,0), limits = c(850,-10)) +
#     scale_y_continuous(name = expression(paste(R~"  [",mg~C~m^{-3}~d^{-1},"] ")), expand = c(0.02,0.02), limits = c(-.5, 10)) +
#     geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
#     coord_flip() +
#     theme_bw() +
#   scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
#                     labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
#   scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
#                       labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
# 
#   
#   Att_plots <- ggplot(data = slopes_final, aes(x = depth, y = resp, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
#     # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
#     geom_point() +
#     geom_line() +
#     # geom_point(shape=21, fill = colors, size=2) +
#     # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
#     # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
#     scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
#     scale_y_continuous(name = expression(paste(Carbon~flux,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,1600)) +
#     coord_flip() +
#     # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
#     theme_bw() +
#     theme(axis.text = element_text(size = 11, colour = "black"),
#           axis.title = element_text(size = 11, colour = "black"),
#           axis.text.x=element_text(angle=0),
#           plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
#     scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
#                       labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
#     scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
#                         labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
#   
#   Att_plots_annual <- ggplot(data = slopes_final, aes(x = depth, y = ANCP, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
#     # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
#     geom_point() +
#     geom_line() +
#     # geom_point(shape=21, fill = colors, size=2) +
#     # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
#     # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
#     scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
#     scale_y_continuous(name = expression(paste(Annual~carbon~flux,"  [",g~C~m^{-2}~y^{-1},"]")), expand = c(0,0),limits = c(-5,250)) +
#     coord_flip() +
#     # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
#     theme_bw() +
#     theme(axis.text = element_text(size = 11, colour = "black"),
#           axis.title = element_text(size = 11, colour = "black"),
#           axis.text.x=element_text(angle=0),
#           plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
#     scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
#                       labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
#     scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
#                         labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
#   
#   
#   layout <- c(
#     area(1,1),
#     area(1,2),
#     area(1,3),
#     area(2,1),
#     area(2,2)
#   )
#   I_plots[[1]] <- I_plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
#                                        axis.title.y.right = element_blank(),axis.text.y.right = element_blank())
#   I_plots[[2]] <- I_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank(),
#                                        axis.text.x = element_blank(), axis.text.y.left = element_blank(),
#                                        axis.title.y.right = element_blank(),axis.text.y.right = element_blank())
#   I_plots[[3]] <- I_plots[[3]] + theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank(),
#                                        axis.title.x = element_blank())
#   I_plots[[4]] <- I_plots[[4]] + theme(axis.title.y.right = element_blank(), axis.text.y.right = element_blank(),
#                                        axis.title.x = element_blank())
#   I_plots[[5]] <- I_plots[[5]] + theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank())
#   
#   finalPlot_I <- wrap_plots(I_plots,design=layout) #+ plot_annotation(tag_levels = "1")
#   ggsave(plot = finalPlot_I,"figures/I_plots_AOU.png", width = 11.5, height = 5, dpi = 300)
#   
#   
#   ggsave(plot = Att_plots,"figures/Att_plots_AOU.png", width = 6, height = 6, dpi = 300)
#   write_csv(slopes_final, "data/slopes_AOU.csv")
#   
#   ggsave(plot = Att_plots_annual,"figures/Att_plots_annual_AOU.png", width = 6, height = 6, dpi = 300)
#   
# 
# }
# 
# 
# AOU_flux(data_DOXY = read_csv("data/DOXY/plotData.csv") %>% arrange(cluster),
#          zp_data = read_csv("data/zp_data_all.csv") %>% arrange(cluster),
#          mld_data = read_csv("data/mld_data_all.csv") %>% arrange(cluster),
#          SIG_data = read_csv("data/Chla_BBP_data/plotData_env.csv") %>% filter(PARAM == "SIG") %>% arrange(cluster),
#          # metrics = read_csv("data/metrics_cluster.csv") %>% arrange(cluster),
#          window = 30)

DOXY_climato_flux(data_DOXY = read_csv("data/DOXY/plotData_climato_DOXY.csv") %>% arrange(cluster),
         zp_data = read_csv("data/zp_data_all.csv") %>% arrange(cluster),
         mld_data = read_csv("data/mld_data_all.csv") %>% arrange(cluster),
         SIG_data = read_csv("data/Chla_BBP_data/plotData_env.csv") %>% filter(PARAM == "SIG") %>% arrange(cluster),
         # metrics = read_csv("data/metrics_cluster.csv") %>% arrange(cluster),
         window = 30)


DOXY_climato_flux <- function(data_DOXY, zp_data, mld_data, SIG_data, window, AOU_sd){
  
  names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
  
  for (i in 1:length(unique(data_DOXY$cluster))){
    
    print(paste("cluster",i))
    
    DOXY_profile_data <- data_DOXY %>%
      filter(cluster == i) %>%
      left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(SIG_data %>% dplyr::select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) 
    # mutate(depths_below_zp = NA)
    # DOXY_profile_data$zp <- round(DOXY_profile_data$zp)
    
    
    depths = DOXY_profile_data$PRES
    dates = DOXY_profile_data$jDay
    month = DOXY_profile_data$month
    zp = DOXY_profile_data$zp
    mld = DOXY_profile_data$mld
    depth_bin = 100
    SIG <- DOXY_profile_data$SIG +1000
    values <- DOXY_profile_data$value * SIG
    depths_below_zp = round(depths-zp)
    
    # definir prof max sous zp pour le cluster
    prof_max <- 1000-max(zp)
    
    if(prof_max>850){
      # analyse
      depth_seq <- c("whole_column", "productive_layer", seq(0, 850-depth_bin, by = depth_bin))
      depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 850, sep = "_"), depth_seq)
      
      # depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, (as.numeric(depth_seq)+depth_bin), sep = "_"), depth_seq)
      
    }else if(prof_max>800){
      # analyse
      depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
      depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 800, sep = "_"), depth_seq)
      # depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, (as.numeric(depth_seq)+depth_bin), sep = "_"), depth_seq)
      
    }else if(prof_max<800){
      # analyse
      depth_seq <- c("whole_column", "productive_layer", seq(0, 650-depth_bin, by = depth_bin))
      depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 650, sep = "_"), depth_seq)}
    # depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, (as.numeric(depth_seq)+depth_bin), sep = "_"), depth_seq)}
    
    
    data_list_O2 <- data.frame(values, depths, dates, zp, mld, depths_below_zp) %>% 
      arrange(depths_below_zp) %>%
      # filter(between(depths_below_zp,0,as.numeric(last(depth_seq))+50)) %>%
      split(f = .$depths_below_zp)
    
    
    # following Su et al 2022
    ZP=round(max(zp))+1
    maxi=as.numeric(last(depth_seq))+50+ZP
    resp_each_depth <- data_list_O2[ZP:maxi] %>%
      llply(function(x){
        metrics <- data.frame(x) %>%
          dplyr::mutate(smoothed_values = smth(values, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>% # smoothed data
          dplyr::summarise(values=values,
                           smoothed_values=smoothed_values,
                           dates=dates,
                           min_y = min(smoothed_values, na.rm=T), # minimum on smoothed data
                           max_y = max(smoothed_values, na.rm=T), # maximum on smoothed data
                           min_x = match(min(smoothed_values, na.rm=T),smoothed_values) %>% as.numeric(), # minimum on smoothed data
                           max_x = match(max(smoothed_values, na.rm=T),smoothed_values) %>% as.numeric(), # maximum on smoothed data
                           Dprod = if(min_x<max_x){Dprod = 365 - max_x + min_x} # productive period between min and max
                           else if(min_x>max_x){Dprod = min_x-max_x},
                           min_supp_max = if(min_x>max_x){min_supp_max = TRUE}
                           else if(min_x<max_x){min_supp_max = FALSE}) %>%
          dplyr::mutate(dates_recal = ifelse(min_supp_max==F & between(dates,0,min_x),dates+365,dates), # recal data if max is before min pretending continuous TS
                        min_x_recal = ifelse(min_supp_max==F,min_x+365,min_x)) %>%
          dplyr::filter(dates_recal<min_x_recal & dates_recal>max_x) %>%
          dplyr::mutate(slopes = summary(lm(values~ dates_recal, data = .))$coefficients[2]* -(1/1.45) *10^-3 *12, # linear regression coeff on unsmoothed data and µmol/m3 O2 to mg C/m3
                        intercept = summary(lm(values~ dates_recal, data = .))$coefficients[1],
                        slope_annual = slopes * Dprod)
      })
    
    # timing de la periode productive
    timings <- resp_each_depth %>% ldply(function(x){
      values_to_plot <- data.frame(x) %>%
        dplyr::summarise(max_x = unique(max_x),
                         min_x = unique(min_x_recal))
    })
    names(timings)[1] <- "PRES"
    timings$PRES <- as.numeric(timings$PRES)
    timings <- timings %>% dplyr::mutate(cluster=i)
    
    # recap des slopes journalières et annuelles
    respiration <- resp_each_depth %>% ldply(function(x){
      values_to_plot <- data.frame(x) %>%
        dplyr::summarise(R = unique(slopes),
                         annual_R = unique(slope_annual))
    })
    names(respiration)[1] <- "depths"
    respiration$depths <- as.numeric(respiration$depths)
    respiration <- respiration %>% dplyr::mutate(cluster=i)
    
    
    ANCP <- data.frame(ANCP=rep(NA,length(depth_layers)), 
                       depth=rep(NA,length(depth_layers)),
                       resp=rep(NA,length(depth_layers)),
                       cluster=rep(i,length(depth_layers)))
    for(l in 1:length(depth_layers)){
      if (grepl("depth", depth_layers[l])) {
        min_depth_to_keep <- str_split(gsub("depth_", "", depth_layers[l]), "_", simplify = T)[,1] %>% as.numeric()
        max_depth_to_keep <- str_split(gsub("depth_", "", depth_layers[l]), "_", simplify = T)[,2] %>% as.numeric()
        index <- which(between(respiration$depths, min_depth_to_keep, max_depth_to_keep))
        ANCP$resp[l] <- (respiration$R[index] * {c( diff(respiration$depths[index]), diff(respiration$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T) #mgC/m2/jour
        ANCP$depth[l] <- min_depth_to_keep
        ANCP$ANCP[l] <- (respiration$annual_R[index] * {c( diff(respiration$depths[index]), diff(respiration$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T)*10^-3 #gC/m2/an
      }
      
      ANCP$Teff <- ANCP$ANCP / ANCP$ANCP[3]
    }
    
    # agreger les slopes ensemble a chaque iterration pour avoir un fichier final
    # avec tout
    
    if(i==1){
      slopes_final <- ANCP
      resp_final <- respiration
      timings_all <- timings
    }
    
    if(i!=1){
      slopes_final <- bind_rows(slopes_final, ANCP)
      resp_final <- bind_rows(resp_final, respiration)
      timings_all <- bind_rows(timings_all, timings)
    }
    
  }
  slopes_final <- drop_na(slopes_final)
  timings_all <- drop_na(timings_all)
  
  resp_each_depth_plot <- ggplot(data = resp_final, aes(x = depths, y = R, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    geom_line() +
    # geom_errorbar(aes(ymin = estimate-sd, ymax = estimate+sd), width = 15) +
    scale_x_reverse(name = "Pressure under Zp [dbar]", expand = c(0,0), limits = c(850,-10)) +
    scale_y_continuous(name = expression(paste(R~"  [",mg~C~m^{-3}~d^{-1},"] ")), expand = c(0.02,0.02), limits = c(-.5, 2.5)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  resp_each_depth_plot_annual <- ggplot(data = resp_final, aes(x = depths, y = annual_R, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    geom_line() +
    # geom_errorbar(aes(ymin = estimate-sd, ymax = estimate+sd), width = 15) +
    scale_x_reverse(name = "Pressure under Zp [dbar]", expand = c(0,0), limits = c(850,-10)) +
    scale_y_continuous(name = expression(paste(R~"  [",mg~C~m^{-3}~d^{-1},"] ")), expand = c(0.02,0.02), limits = c(-25, 300)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  
  Att_plots <- ggplot(data = slopes_final, aes(x = depth, y = resp, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
    geom_point() +
    geom_line() +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
    # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
    scale_y_continuous(name = expression(paste(Carbon~flux,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,450)) +
    coord_flip() +
    # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
    theme_bw() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.text.x=element_text(angle=0),
          plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  Att_plots_annual <- ggplot(data = slopes_final, aes(x = depth, y = ANCP, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
    geom_point() +
    geom_line() +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
    # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
    scale_y_continuous(name = expression(paste(Annual~carbon~flux,"  [",g~C~m^{-2}~y^{-1},"]")), expand = c(0,0),limits = c(-5,80)) +
    coord_flip() +
    # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
    theme_bw() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.text.x=element_text(angle=0),
          plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  Teff_plots_annual <- ggplot(data = slopes_final, aes(x = depth, y = Teff, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
    geom_point() +
    geom_line() +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
    # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
    scale_y_continuous(name = expression(paste(Teff,"  [%]")), expand = c(0,0),limits = c(0,1)) +
    coord_flip() +
    # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
    theme_bw() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.text.x=element_text(angle=0),
          plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))

  
  ggsave(plot = Att_plots,"figures/Att_plots_DOXY_climato.png", width = 6, height = 6, dpi = 300)
  write_csv(slopes_final, "data/slopes/slopes_DOXY_climato.csv")
  
  ggsave(plot = Att_plots_annual,"figures/Att_plots_annual_DOXY_climato.png", width = 6, height = 6, dpi = 300)
  
  
}




O2_POC <- function(slopes_DOXY = read_csv("data/slopes_DOXY_climato.csv", show_col_types = F) %>% arrange(depth),
                   slopes_POC = read_csv("data/slopes_POC.csv", show_col_types = F) %>% arrange(DEPTH)){
  
  Rsub_Att <- data.frame(slopes_POC = slopes_POC$slope,
               # slope_min_POC = slopes_POC$slope_min,
               # slope_max_POC = slopes_POC$slope_max,
               slopes_AOU = slopes_DOXY$resp,
               # slope_min_AOU = slopes_DOXY$slope_min,
               # slope_max_AOU = slopes_DOXY$slope_max,
               cluster = slopes_DOXY$cluster,
               prof = slopes_POC$DEPTH,
               slope_annual_AOU = slopes_DOXY$ANCP,
               slope_annual_POC = slopes_POC$slope_annual * 10^-3)
  
  Att_plots_annual_DOXY <- ggplot(data = Rsub_Att, aes(x = prof, y = slope_annual_AOU, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
    geom_point() +
    geom_line() +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
    # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
    scale_y_continuous(name = expression(paste(Annual~carbon~flux~TOT,"  [",g~C~m^{-2}~y^{-1},"]")), expand = c(0,0),limits = c(-5,80)) +
    coord_flip() +
    # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
    theme_bw() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.text.x=element_text(angle=0),
          plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  Att_plots_annual_POC <- ggplot(data = Rsub_Att, aes(x = prof, y = slope_annual_POC, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
    # geom_ribbon(aes(ymin = slope_min, ymax = slope_max), alpha = 0.3, linetype = 0) +
    geom_point() +
    geom_line() +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
    # geom_line(data = regression_data, aes(x = z, y = slope), size = 0.5) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(850,0)) +
    scale_y_continuous(name = expression(paste(Annual~carbon~flux~SMALL,"  [",g~C~m^{-2}~y^{-1},"]")), expand = c(0,0),limits = c(-.5,10)) +
    coord_flip() +
    # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
    theme_bw() +
    theme(axis.text = element_text(size = 11, colour = "black"),
          axis.title = element_text(size = 11, colour = "black"),
          axis.text.x=element_text(angle=0),
          plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
    scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
    scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                        labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))
  
  Rsub_Att <- Rsub_Att%>%
    mutate(ratio = slope_annual_POC/slope_annual_AOU * 100) %>% arrange(cluster) 
    # mutate(clust_name = c(rep("(1) Austral HCB",8),rep("(2) Austral MCB",8),rep("(3) STZ",8),rep("(4) PN",8),rep("(5) AN",8))) %>%

  
  ratio_plot_annual <- ggplot(data = Rsub_Att, aes(x = as.numeric(prof), y = ratio, color = factor(cluster))) +
    geom_point(aes(colour = factor(cluster))) +
    geom_line(aes(colour = factor(cluster))) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0.01,0.01), limits = c(750,0)) +
    scale_y_continuous(name = paste("contribution des petites particules à l'export TOT annuel [%]"), expand = c(0.01,0.01),limits = c(0,55)) +
    coord_flip() +
    theme_bw() +
    scale_color_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                       labels = c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")) +
    theme(#axis.title.x = element_blank(),
      #axis.text.x = element_blank(),
      legend.title = element_blank())
  
  ggsave(plot = ratio_plot_annual,"figures/small_TOT_VS_depths_annual.png", width = 6, height = 6, dpi = 300)
}

O2_POC()

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


# Check that it doesn't match any non-number
numbers_only <- function(x) !grepl("\\D", x)
