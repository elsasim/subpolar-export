
library(factoextra)
library(tidyverse)
library(lubridate)
library(zoo)
library(RColorBrewer)
library(cowplot)
# library(dplyr)
# library(tidyr)
# library(readr)
# library(stringr)
library(data.table)
library(patchwork)

library(plyr)
#library(dplyr)
#library(ggplot2)
#library(readr)
#library(stringr)
library(smoother)
#library(tidyr)
library(scales)
library(grid)
library(numDeriv)
library(stats)

library(broom)

plotData_DOXY <- read_csv("data/all/DOXY_data_all.csv") 
plotData_MLD <- read_csv("data/all/mld_data_all.csv")
plotData_SIG <- read_csv("data/all/SIG_data_all.csv")

mask_data <- read_csv("data/all/mask_data_all.csv") %>% arrange(cluster)


# observer la diminution d'O2 à une profondeur constante et en déduire des taux de respiration sur la profondeur puis une respiration intégrée
#####################################
Hennon_2016 <-function(plotData_DOXY, plotData_MLD, mask_data, window = 10){
  
  O2_plots <- list()
  R_plots <- list()
  regression_plot_70 <- list()
  regression_plot_150 <- list()
  regression_plot_250 <- list()
  
  Rsub <- vector()
  chla <- vector()
  POC <- vector()
  max_POC <- vector()
  max_chla <- vector()
  
  for(i in 1:length(unique(plotData_DOXY$cluster))){
    
      DOXY_profile_data <- plotData_DOXY %>%
        filter(PARAM == "DOXY", cluster == i) %>%
        left_join(plotData_MLD %>% select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
        left_join(plotData_SIG %>% select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster"))
      
      rates <- data.frame(PRES = seq(70, 1000, by=10), dO2 = NA)
    
      for (j in 7:100){
          data_isodepth <- filter(DOXY_profile_data, PRES == (j*10)) %>%
            filter(PRES > mld) %>%
            mutate(smth_values = smth(value, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>%
            select(jDay, value, smth_values)
          
          regression <- lm(data = data_isodepth, formula = value ~ jDay)
          rates$dO2[j-6] <- regression$coefficients[[2]]
          
         if(j==7){
            regression_plot_70[[i]] <- ggplot(data=data_isodepth, aes(x=jDay, y=smth_values)) +
              geom_point(size = 1, color = "black") +
              scale_x_continuous(name = "relative Julian Date", expand = c(0,0), limits = c(0,365)) +
              scale_y_continuous(name = expression(paste("[",O[2],"]"~µmol~O[2]~kg^{-1})), expand = c(0,0), limits = c(240,310)) +
              theme_bw() +
              geom_abline(slope = regression$coefficients[[2]], intercept = regression$coefficients[[1]], color = "grey")
         }
          if(j==15){
            regression_plot_150[[i]] <- ggplot(data=data_isodepth, aes(x=jDay, y=smth_values)) +
              geom_point(size = 1, color = "black") +
              scale_x_continuous(name = "relative Julian Date", expand = c(0,0), limits = c(0,365)) +
              scale_y_continuous(name = expression(paste("[",O[2],"]"~µmol~O[2]~kg^{-1})), expand = c(0,0), limits = c(125,310)) +
              theme_bw() +
              geom_abline(slope = regression$coefficients[[2]], intercept = regression$coefficients[[1]], color = "grey")
          }
          
          if(j==25){
            regression_plot_250[[i]] <- ggplot(data=data_isodepth, aes(x=jDay, y=smth_values)) +
              geom_point(size = 1, color = "black") +
              scale_x_continuous(name = "relative Julian Date", expand = c(0,0), limits = c(0,365)) +
              scale_y_continuous(name = expression(paste("[",O[2],"]"~µmol~O[2]~kg^{-1})), expand = c(0,0), limits = c(60,310)) +
              theme_bw() +
              geom_abline(slope = regression$coefficients[[2]], intercept = regression$coefficients[[1]], color = "grey")
          }
          
      }
      sigma <- mean(DOXY_profile_data$SIG) + 1000    
      rates$R <- -(1/1.45)*rates$dO2 * 365 *10^-3 *12 *sigma  #O2 in C, days in years and µmol in mg and kg-1 in m-3
  
      O2_plots[[i]] <- ggplot(data = rates, aes(x = PRES, y = dO2)) +
        geom_line(size = 1, color = "black") +
        scale_x_reverse(name = "Pressure [dbar]", expand = c(0,0), limits = c(1000,0)) +
        scale_y_continuous(name = expression(paste(dO[2]/dt,"  [",µmol~O[2]~kg^{-1}~d^{-1},"]")), expand = c(0,0),limits = c(-0.25,0.25)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
        coord_flip() +
        theme_bw()
      
      
      R_plots[[i]] <- ggplot(data = rates, aes(x = PRES, y = R)) +
        geom_line(size = 1, color = "black") +
        scale_x_reverse(name = "Pressure [dbar]", expand = c(0,0), limits = c(1000,0)) +
        scale_y_continuous(name = expression(paste(net~respiration~"  [",mg~C~m^{-3}~y^{-1},"] ")), expand = c(0.1,0.1)) +
        geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
        coord_flip() +
        theme_bw()
      
      
      ################
      # Rsub = depth integrated respiration
      
      resp <- rates$R  
      PRES <- rates$PRES

      Rsub[i] <- abs(PRES[2]-PRES[1])*abs(resp[2]-resp[1])
      k=2
      
      while(resp[k] > 0 & k<=93 ){
        Rsub[i] <- Rsub[i] + abs(PRES[k+1]-PRES[k])*abs(resp[k+1]-resp[k])
        k=k+1
      }
      
     # Rsub[i] <- Rsub[i]#convert mmol in mol/m-2/y-1
       
     ######## CHL mask
    CHL_data <- filter(mask_data, PARAM=="CHL_masked" & cluster==i)
     dates <- CHL_data$jDay
     CHL_mask <- CHL_data$value
     
     chla[i] <- mean(CHL_mask)
     max_chla[i] <- max(CHL_mask)
     
     ### POC mask
     POC_data <- filter(mask_data, PARAM=="POC_masked" & cluster==i)
     values_iPOC_mask <- vector() %>% as.numeric()
     dates <- POC_data$jDay
     POC_mask <- POC_data$value #in mgC
     
     # for (l in 1:length(unique(dates))){
     #   range <- which(dates==unique(dates)[l])
     #   values_iPOC_mask[l] <- sum(POC_mask[range]*10^-3)
     # }
     POC[i] <- mean(POC_mask)
     max_POC[i] <- max(POC_mask)
     
  }
  
  comparision <- data.frame(chla = chla, max_chla = max_chla, POC = POC, max_POC = max_POC, Rsub = Rsub) #data frame with all data
  colors = c("#CC0000","#FF6666","#FFCC00","#33CCFF", "#0066CC")
  
  chla_plot <- ggplot(data = comparision, aes(x = chla, y = Rsub)) +
    geom_point(shape=21, fill = colors, size = 3) +
    scale_x_continuous(name = expression(paste("Chla  [",mg~m^{-3},"]")), expand = c(0.01,0.01)) +
    scale_y_continuous(name = expression(paste(R[sub]~~"[",mg~C~m^{-2}~y^{-1},"]")), expand = c(0.1,0.1), limits = c(0.1,0.6)) +
    theme_bw()
    # labels = c("Austral_HCB (1)","Austral_MCB (2)","Austral_LCB (3)","PN (4)","AN (5)"))
  
  max_chla_plot <- ggplot(data = comparision, aes(x = max_chla, y = Rsub)) +
    geom_point(shape=21, fill = colors, size = 3) +
    scale_x_continuous(name = expression(paste("Max Chla  [",mg~m^{-3},"]")), expand = c(0.01,0.01)) +
    scale_y_continuous(name = expression(paste(R[sub]~~"[",mg~C~m^{-2}~y^{-1},"]")), expand = c(0.1,0.1), limits = c(0.1,0.6)) +
    theme_bw()
    # scale_fill_manual(values = c("#CC0000","#FF6666","#FFCC00","#33CCFF", "#0066CC"),
    #                   labels = c("Austral_HCB (1)","Austral_MCB (2)","Austral_LCB (3)","PN (4)","AN (5)"))
  
  POC_plot <- ggplot(data = comparision, aes(x = POC, y = Rsub)) +
    geom_point(shape=21, fill = colors, size = 3) +
    scale_x_continuous(name = expression(paste(iPOC,"  [",mg~C~m^{-2},"]")), expand = c(0.1,0.1)) +
    scale_y_continuous(name = expression(paste(R[sub]~~"[",mg~C~m^{-2}~y^{-1},"]")), expand = c(0.1,0.1), limits = c(0.1,0.6)) +
    theme_bw()
  # labels = c("Austral_HCB (1)","Austral_MCB (2)","Austral_LCB (3)","PN (4)","AN (5)"))
  
  max_POC_plot <- ggplot(data = comparision, aes(x = max_POC, y = Rsub)) +
    geom_point(shape=21, fill = colors, size = 3) +
    scale_x_continuous(name = expression(paste(Max~iPOC,"  [",mg~C~m^{-2},"]")), expand = c(0.1,0.1)) +
    scale_y_continuous(name = expression(paste(R[sub]~~"[",mg~C~m^{-2}~y^{-1},"]")), expand = c(0.1,0.1), limits = c(0.1,0.6)) +
    theme_bw()
  
  CHL_POC <- ggplot(data = comparision, aes(x = POC, y = chla)) +
    geom_point(shape=21, fill = colors, size = 3) +
    scale_x_continuous(name = expression(paste(iPOC,"  [",mg~C~m^{-2},"]")), expand = c(0.1,0.1)) +
    scale_y_continuous(name = expression(paste("Chla  [",mg~m^{-3},"]")), expand = c(0.1,0.1)) +
    theme_bw()
  
  Max_CHL_POC <- ggplot(data = comparision, aes(x = max_POC, y = max_chla)) +
    geom_point(shape=21, fill = colors, size = 3) +
    scale_x_continuous(name = expression(paste(Max~iPOC,"  [",mg~C~m^{-2},"]")), expand = c(0.1,0.1)) +
    scale_y_continuous(name = expression(paste("Max Chla  [",mg~m^{-3},"]")), expand = c(0.1,0.1)) +
    theme_bw()
  
  O2_plots[[1]] <- O2_plots[[1]] + theme(axis.title.x = element_blank())
  O2_plots[[2]] <- O2_plots[[2]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  O2_plots[[3]] <- O2_plots[[3]] + theme(axis.title.y = element_blank())
  #O2_plots[[4]] <- O2_plots[[4]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  O2_plots[[5]] <- O2_plots[[5]] + theme(axis.title.y = element_blank())
  
  
  R_plots[[1]] <- R_plots[[1]] + theme(axis.title.x = element_blank())
  R_plots[[2]] <- R_plots[[2]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  R_plots[[3]] <- R_plots[[3]] + theme(axis.title.y = element_blank())
  #R_plots[[4]] <- R_plots[[4]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  R_plots[[5]] <- R_plots[[5]] + theme(axis.title.y = element_blank())
  
  final_plot_70 <- regression_plot_70[[1]] + regression_plot_70[[2]] + regression_plot_70[[3]]+ regression_plot_70[[4]]+ regression_plot_70[[5]]
  final_plot_150 <- regression_plot_150[[1]] + regression_plot_150[[2]] + regression_plot_150[[3]]+ regression_plot_150[[4]]+ regression_plot_150[[5]]
  final_plot_250 <- regression_plot_250[[1]] + regression_plot_250[[2]] + regression_plot_250[[3]]+ regression_plot_250[[4]]+ regression_plot_250[[5]]
  
  layout <- c(
    area(1,1),
    area(1,2),
    area(1,3),
    area(2,1),
    area(2,2)
  )
  
  finalPlot_O2 <- wrap_plots(O2_plots,design=layout) 
  # ggsave(plot = finalPlot_O2,"figures/O2_plots.png", width = 11, height = 8, dpi = 300)
  
  finalPlot_R <- wrap_plots(R_plots,design=layout) 
  # ggsave(plot = finalPlot_R,"figures/R_plots.png", width = 11, height = 8, dpi = 300)
  
  # ggsave(plot = final_plot_70,"figures/O2_iso70_plots.png", width = 11, height = 8, dpi = 300)
  # ggsave(plot = final_plot_150,"figures/O2_iso150_plots.png", width = 11, height = 8, dpi = 300)
  # ggsave(plot = final_plot_250,"figures/O2_iso250_plots.png", width = 11, height = 8, dpi = 300)
  
  ggsave(plot = max_chla_plot,"figures/Rsub_max_chla_plot.png", width = 5, height = 3.6, dpi = 300)
  ggsave(plot = chla_plot,"figures/Rsub_chla_plot.png", width = 5, height = 3.6, dpi = 300)
  
  ggsave(plot = max_POC_plot,"figures/Rsub_max_POC_plot.png", width = 5, height = 3.6, dpi = 300)
  ggsave(plot = POC_plot,"figures/Rsub_POC_plot.png", width = 5, height = 3.6, dpi = 300)
  
  ggsave(plot = Max_CHL_POC,"figures/Max_CHL_POC_plot.png", width = 5, height = 3.6, dpi = 300)
  ggsave(plot = CHL_POC,"figures/CHL_POC_plot.png", width = 5, height = 3.6, dpi = 300)
}

Hennon_2016(plotData_DOXY, plotData_MLD, mask_data)

##########################################################################

plotData_DOXY <- read_csv("data/all/DOXY_data_all.csv") 
plotData_MLD <- read_csv("data/all/mld_data_all.csv")
plotData_SIG <- read_csv("data/all/SIG_data_all.csv")
zp_data <- read_csv("data/all/zp_data_all.csv") %>% arrange(cluster)
mld_data <- read_csv("data/all/mld_data_all.csv") %>% arrange(cluster)

mask_data <- read_csv("data/all/mask_data_all.csv") %>% arrange(cluster)

Hennon_V2 <- function(plotData_DOXY, plotData_MLD, mask_data, window=30){
  
  R_plots <- list()
  
  for (i in 1:5){
    
    print(paste("cluster", i))
    
    DOXY_profile_data <- plotData_DOXY %>%
      filter(PARAM == "DOXY", cluster == i) %>%
      left_join(zp_data %>% select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(mld_data %>% select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(plotData_SIG %>% select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) 
    
    
    POC_data <- filter(mask_data, PARAM=="POC_masked" & cluster==i)
    CHL_data <- filter(mask_data, PARAM=="CHL_masked" & cluster==i)
    
    values <- DOXY_profile_data$value
    
    depths = DOXY_profile_data$PRES
    dates = DOXY_profile_data$jDay
    zp = DOXY_profile_data$zp
    mld = DOXY_profile_data$mld
    depth_bin = 150
    
    depth_seq <- c("whole_column", "productive_layer", seq(0, 1000-depth_bin, by = depth_bin))
    # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
    # 
    depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 1000, sep = "_"), depth_seq)
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
          index <- which(between(df$depths_below_zp,  min_depth_to_keep, max_depth_to_keep))
        }
        #(df$values[index] * {c( diff(df$depths[index]), diff(df$depths[index]) %>% dplyr::last() )}) %>% sum(na.rm = T) # The integration
        mean(df$values[index])
      })

      data.frame(df[1,], layer_name = depth_layers, DEPTH = depth_seq, I = integrated_values, row.names = NULL)
      
    }, .parallel = FALSE, .id = NULL)
    
    values_per_layers <- values_per_layers %>%
      dplyr::group_by(layer_name) %>%
      dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T)) %>%
      ungroup()
    
    depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
    colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")

    if(i==1){
      max_mld <- 189
    }
    if(i==2){
      max_mld <- 178
    }
    if(i==3){
      max_mld <- 117
    }
    if(i==4){
      max_mld <- 89
    }
    if(i==5){
      max_mld <- 149
    }
    
    values_per_layers$DEPTH <- as.numeric(values_per_layers$DEPTH)
    
    my_slope <- values_per_layers %>%
      filter(DEPTH == "0"|DEPTH == "150"|DEPTH == "300"|DEPTH == "450"|DEPTH == "600"|DEPTH == "750") %>%
      filter(dates >= max_mld) %>%
      select(dates = dates, O2 = smoothed_I, DEPTH = DEPTH) %>%
      group_by(DEPTH) %>% do(model = lm(O2~ dates, data = .)) %>% ungroup() %>% #groupe les données par "lvl" et applique un lm pour chaque "lvl". La fonction sors une liste pour chaque "lvl"
      transmute(DEPTH, coef = map(model, tidy)) %>% #la fonction tidy de broom permet de récupérer les infos de la list comme le coef std.error etc... Plour le R² modifier par la fonction "glance"
      unnest(coef) %>%
      filter(term == "dates") %>% #la fonction précédente retourne tout en une seule colonne. Unnest la sépare en plusieur colonne
      #optionel : filtrer les lignes relatives à l'intercept pour ne travailler que sur Y
      select(DEPTH,term, estimate, std.error) #optionnel: sélectionner que les "lvl" et intercept associés
  
    sigma <- mean(DOXY_profile_data$SIG) + 1000
    
    my_slope$estimate <- -(1/1.45)*my_slope$estimate *10^-3 *12 *sigma

   my_slope <-  my_slope %>% 
      mutate(sd = (1/1.45)*std.error *10^-3 *12 *sigma)
    
    R_plots[[i]] <- ggplot(data = my_slope, aes(x = DEPTH, y = estimate)) +
      geom_point(shape=21, fill = colors, size=2) +
      geom_errorbar(aes(ymin = estimate-sd, ymax = estimate+sd), width = 15) +
      scale_x_reverse(name = "Pressure under Zp [dbar]", expand = c(0,0), limits = c(1000,-10)) +
      scale_y_continuous(name = expression(paste(net~respiration~"  [",mg~C~m^{-3}~d^{-1},"] ")), expand = c(0.02,0.02), limits = c(-0.15, 0.45)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
      coord_flip() +
      theme_bw()
    
   slope <- my_slope$estimate
   DEPTH <- my_slope$DEPTH
   Rsub <- vector()
   for (j in 1:6){
     Rsub[j] <- (DEPTH[j+1] - DEPTH[j])*(slope[j]-slope[j+1])
   }
    
  }
  
  layout <- c(
    area(1,1),
    area(1,2),
    area(1,3),
    area(2,1),
    area(2,2)
  )
  
  R_plots[[1]] <- R_plots[[1]] + theme(axis.title.x = element_blank())
  R_plots[[2]] <- R_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank())
  R_plots[[3]] <- R_plots[[3]] + theme(axis.title.y.left = element_blank())
  #R_plots[[4]] <- R_plots[[4]] 
  R_plots[[5]] <- R_plots[[5]] + theme(axis.title.y.left = element_blank())
  
  finalPlot_R <- wrap_plots(R_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  
  ggsave(plot = finalPlot_R,"figures/R_plot.png",  width = 11, height = 8, dpi = 300)
  
}
Hennon_V2(plotData_DOXY, plotData_MLD, mask_data, window=30)
  