




#### Load packages
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
library(raster)
library(tidyverse)



#################################################################################

Dall_Olmo_2014_processing_BBP <- function(profile_data, zp_data, mld_data, mask_data, POC_data, window=30) {
  
  # plot_theme <-     theme(text = element_text(size=25, colour = "black"),
  #                         plot.title = element_text(hjust = 0.5),
  #                         panel.grid.major = element_blank(),
  #                         panel.grid.minor = element_blank(),
  #                         panel.background = element_blank(),
  #                         panel.border = element_rect(linetype = "solid", fill = NA),
  #                         axis.text = element_text(size = 25, colour = "black"),
  #                         axis.title = element_text(size = 25, colour = "black"),
  #                         axis.text.x=element_text(angle=0))
  
  I_plots <- list()
  Ez_plots <- list()
  Rp_plots <- list()
  Att_plots <- list()
  Instant_flux_plots <- list()
  Instant_norm_flux_plots <- list()
  sinking_speed_plots <- list()
  
  R_plots <- list()
  Rsub_plots <- list()
  Rsub_Att <- list()
  Rsub_Att_plots <- list()
  Rsub_0_100 <- vector()
  Rsub_100_200 <- vector()
  Rsub_200_300 <- vector()
  Rsub_300_400 <- vector()
  Rsub_400_500 <- vector()
  Rsub_500_600 <- vector()
  Rsub_600_700 <- vector()
  Rsub_700_800 <- vector()
  
  regression_Rsub_Att <- list()
  max_mld <- vector()
  O2 <- list()
  Martin_exp <- vector()
  names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
  # incertitudes_C <- read_csv("data/all/incertitudes_C.csv")
  # incertitudes_O <- read_csv("data/all/incertitudes_O2.csv")
  # 
  # my_slope <- data.frame(depths_below_zp=rep(NA,954),term=rep(NA,954), estimate=rep(NA,954), std.error =rep(NA,954))
  # my_slope <- data.frame(depths_below_zp=rep(NA,6),term=rep(NA,6), estimate=rep(NA,6), std.error =rep(NA,6))
  
  for (i in 1:length(unique(profile_data$cluster))){
    
    print(paste("cluster", i))
    
    bbp_profile_data <- profile_data %>%
      filter(PARAM == "BBP700", cluster == i) %>%
      left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(mask_data %>% dplyr::select(jDay, PRES, Cphyto_mask = value, cluster), by = c("jDay", "PRES", "cluster"))%>%
      # left_join(POC_data %>% dplyr::select(jDay, PRES, POC = value, cluster), by = c("jDay", "PRES", "cluster")) 
      left_join(POC_data %>% filter( PARAM == "POC_Koestner", cluster == i) %>% dplyr::select(jDay, PRES,POC = value, cluster), by = c("jDay", "PRES", "cluster"))
    
    values <- bbp_profile_data$POC*10^-3 # mg.m-3 en g.m-3
    
    depths = bbp_profile_data$PRES
    dates = bbp_profile_data$jDay
    month = bbp_profile_data$month
    zp = bbp_profile_data$zp
    mld = bbp_profile_data$mld
    Cphyto_mask = bbp_profile_data$Cphyto_mask
    depth_bin = 100
    
    #############################################################################
    ############################ PLOT MASK FOR METRIC ###########################
    #############################################################################
    #############################################################################
    values_iCphyto_mask <- vector() %>% as.numeric()
    for (k in 1:length(unique(dates))){
      range <- which(dates==unique(dates)[k])
      values_iCphyto_mask[k] <- sum(Cphyto_mask[range]*10^-3)
    }
    smth_values_iCphyto_mask <- smth(values_iCphyto_mask, window = window, alpha = 2.5, method = "gaussian", tails = T)
    
    Ez_Cphyto_mask <- c( NA, diff(smth_values_iCphyto_mask) / {diff(unique(dates))})
    rp_values <- Ez_Cphyto_mask/smth_values_iCphyto_mask
    
    apex_x <- match(max(smth_values_iCphyto_mask),smth_values_iCphyto_mask) 
    apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == apex_x)])
    # onset_x <- match(min(smth_values_iPOC_mask[0:apex_x]),smth_values_iPOC_mask[0:apex_x])
    onset_x <- match(min(smth_values_iCphyto_mask),smth_values_iCphyto_mask)
    onset_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == onset_x)])
    
    range_pos <- vector()
    
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
      
      climax_x <- match(max(range_pos[0:apex_x], na.rm = T) , range_pos)
      climax_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == climax_x)])
      
    }
    
    apex_y <- max(smth_values_iCphyto_mask)
    climax_y <- smth_values_iCphyto_mask[climax_x]
    onset_y <- min(smth_values_iCphyto_mask[0:apex_x])
    
    metrics <- c(onset_x, climax_x, apex_x) %>% as.data.frame()
    write_csv(metrics, paste("data/metrics_cluster",i,".csv", sep=""))
    
    ###########
    #####
    scaling_factor_2 <- sec_axis_adjustement_factors(c(0,8), c(-0.02,0.03))
    
    ### accumulation rates
    
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
    
    #############################################################################
    ############################ ATTENUATION PLOT ###############################
    #############################################################################
    #############################################################################
    
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
    

    whole_column <- values_per_layers %>%
      filter(DEPTH == "whole_column")
    
    iPOC_whole_column <- whole_column$smoothed_I
    
    productive_layer <- values_per_layers %>%
      filter(DEPTH == "productive_layer")
    iPOC_productive_layer <- productive_layer$smoothed_I
    
    depth_0_800 <- values_per_layers %>%
      filter(DEPTH == "0")
    iPOC_0_800 <- depth_0_800$smoothed_I
    
    depth_100_800 <- values_per_layers %>%
      filter(DEPTH == "100")
    iPOC_100_800 <- depth_100_800$smoothed_I
    
    depth_200_800 <- values_per_layers %>%
      filter(DEPTH == "200")
    iPOC_200_800 <- depth_200_800$smoothed_I
    
    depth_300_800 <- values_per_layers %>%
      filter(DEPTH == "300")
    iPOC_300_800 <- depth_300_800$smoothed_I
    
    depth_400_800 <- values_per_layers %>%
      filter(DEPTH == "400")
    iPOC_400_800 <- depth_400_800$smoothed_I
    
    depth_500_800 <- values_per_layers %>%
      filter(DEPTH == "500")
    iPOC_500_800 <- depth_500_800$smoothed_I
    
    depth_600_800 <- values_per_layers %>%
      filter(DEPTH == "600")
    iPOC_600_800 <- depth_600_800$smoothed_I
    
    depth_700_800 <- values_per_layers %>%
      filter(DEPTH == "700")
    iPOC_700_800 <- depth_700_800$smoothed_I
    
    ######### productive
    # max_productive_layer_x <- match(max(iPOC_productive_layer),iPOC_productive_layer)
    # max_productive_layer_y <- iPOC_productive_layer[max_productive_layer_x]
    # min_productive_layer_x <- match(min(iPOC_productive_layer[0: max_productive_layer_x]),iPOC_productive_layer[0: max_productive_layer_x])
    # min_productive_layer_y <- iPOC_productive_layer[min_productive_layer_x]
    
    # slope_productive_layer <- (max_productive_layer_y - min_productive_layer_y) / (max_productive_layer_x - min_productive_layer_x) * (max_productive_layer_x-min_productive_layer_x) *10^3
    
    ######### 0
    min_0_800_x <- match(min(iPOC_0_800[0:apex_x]),iPOC_0_800[0:apex_x])
    min_0_800_y <- iPOC_0_800[min_0_800_x]
    max_0_800_x <- match(max(iPOC_0_800),iPOC_0_800)
    max_0_800_y <- iPOC_0_800[max_0_800_x]
    slope_0_800 <- (max_0_800_y - min_0_800_y) / (max_0_800_x - min_0_800_x) 
    # slope_0_800_export <- (max_0_800_y - min_0_800_y) / (max_0_800_x - min_0_800_x) * (max_0_800_x-min_0_800_x) *10^3
    
    ######### 150
    min_100_800_x <- match(min(iPOC_100_800[0:apex_x]),iPOC_100_800[0:apex_x])
    min_100_800_y <- iPOC_100_800[min_100_800_x]
    max_100_800_x <- match(max(iPOC_100_800),iPOC_100_800)
    max_100_800_y <- iPOC_100_800[max_100_800_x]
    slope_100_800 <- (max_100_800_y - min_100_800_y) / (max_100_800_x - min_100_800_x)
    ######### 300
    min_200_800_x <- match(min(iPOC_200_800[0:apex_x]),iPOC_200_800[0:apex_x])
    min_200_800_y <- iPOC_200_800[min_200_800_x]
    max_200_800_x <- match(max(iPOC_200_800),iPOC_200_800)
    max_200_800_y <- iPOC_200_800[max_200_800_x]
    slope_200_800 <- (max_200_800_y - min_200_800_y) / (max_200_800_x - min_200_800_x)
    ######### 450
    min_300_800_x <- match(min(iPOC_300_800[0:apex_x]),iPOC_300_800[0:apex_x])
    min_300_800_y <- iPOC_300_800[min_300_800_x]
    max_300_800_x <- match(max(iPOC_300_800),iPOC_300_800)
    max_300_800_y <- iPOC_300_800[max_300_800_x]
    slope_300_800 <- (max_300_800_y - min_300_800_y) / (max_300_800_x - min_300_800_x)
    ######### 600
    min_400_800_x <- match(min(iPOC_400_800[0:apex_x]),iPOC_400_800[0:apex_x])
    min_400_800_y <- iPOC_400_800[min_400_800_x]
    max_400_800_x <- match(max(iPOC_400_800),iPOC_400_800)
    max_400_800_y <- iPOC_400_800[max_400_800_x]
    slope_400_800 <- (max_400_800_y - min_400_800_y) / (max_400_800_x - min_400_800_x)
    ######### 750
    min_500_800_x <- match(min(iPOC_500_800[0:apex_x]),iPOC_500_800[0:apex_x])
    min_500_800_y <- iPOC_500_800[min_500_800_x]
    max_500_800_x <- match(max(iPOC_500_800),iPOC_500_800)
    max_500_800_y <- iPOC_500_800[max_500_800_x]
    slope_500_800 <- (max_500_800_y - min_500_800_y) / (max_500_800_x - min_500_800_x)
    
    min_600_800_x <- match(min(iPOC_600_800[0:apex_x]),iPOC_600_800[0:apex_x])
    min_600_800_y <- iPOC_600_800[min_600_800_x]
    max_600_800_x <- match(max(iPOC_600_800),iPOC_600_800)
    max_600_800_y <- iPOC_600_800[max_600_800_x]
    slope_600_800 <- (max_600_800_y - min_600_800_y) / (max_600_800_x - min_600_800_x)
    
    min_700_800_x <- match(min(iPOC_700_800[0:apex_x]),iPOC_700_800[0:apex_x])
    min_700_800_y <- iPOC_700_800[min_700_800_x]
    max_700_800_x <- match(max(iPOC_700_800),iPOC_700_800)
    max_700_800_y <- iPOC_700_800[max_700_800_x]
    slope_700_800 <- (max_700_800_y - min_700_800_y) / (max_700_800_x - min_700_800_x)
    
    
    # mettre tout dans un dataframe et passer en mg.m3
    slopes <- data.frame(flux = c(slope_0_800*10^3,slope_100_800*10^3,
                                  slope_200_800*10^3,slope_300_800*10^3,
                                  slope_400_800*10^3,slope_500_800*10^3,
                                  slope_600_800*10^3,slope_700_800*10^3))%>%
      mutate(i_depths = depth_nb+1)

    ###############################################
    incertitudes <- data.frame(cluster=i, layer = depth_layers[3:10], slope_min = rep(NA,8), slope_max = rep(NA,8))
    layer_name <- paste(depth_seq, 800, sep = "_")
    
    for (l in 3:10){
      print(paste("couche", layer_name[l]))
    min_y <- get(paste("min_",layer_name[l],"_y", sep=""))
    max_y <- get(paste("max_",layer_name[l],"_y", sep=""))
    min_x <- get(paste("min_",layer_name[l],"_x", sep=""))
    max_x <- get(paste("max_",layer_name[l],"_x", sep=""))
    sd_tot <- incertitudes_tot_POC(POC_sd, window=30, i)[[l-2]]
    
    if(min_x<max_x){
      # pente minimum
      new_min_y <- min_y + sd_tot[min_x]
      new_max_y <- max_y - sd_tot[max_x]
      slope_min <- (new_max_y - new_min_y) / (max_x - min_x) 
      
      # pente maximum
      new_min_y <- min_y - sd_tot[min_x]
      new_max_y <- max_y + sd_tot[max_x]
      slope_max <- (new_max_y - new_min_y) / (max_x - min_x) 
    }
    
    if(min_x>max_x){
      # pente minimum
      new_min_y <- min_y + sd_tot[min_x]
      new_max_y <- max_y - sd_tot[max_x]
      slope_min <- (new_max_y - new_min_y) / ((365+max_x) - min_x)
      
      # pente maximum
      new_min_y <- min_y - sd_tot[min_x]
      new_max_y <- max_y + sd_tot[max_x]
      slope_max <- (new_max_y - new_min_y) / ((365+max_x) - min_x)
    }
    
    incertitudes$slope_min[l-2] <- slope_min*10^3
    incertitudes$slope_max[l-2] <- slope_max*10^3
    }
    
    
    ###############################################
    slopes <- slopes %>% mutate(slope_min = incertitudes$slope_min,
                                slope_max = incertitudes$slope_max,
                                diff = slope_max - flux)
    
    if(i==5){
      slopes <- slopes[-(7:8),]
      slopes <- slopes %>%
        mutate(prof = c(50,150,250,350,450,550))
      # incert <- filter(incertitudes_C, cluster==i) %>%
      #   mutate(prof = c(50,150,250,350,450,550))
    }
    
    if(i!=5){
      slopes <- slopes %>%
        mutate(prof = c(50,150,250,350,450,550,650,750))
      # incert <- filter(incertitudes_C, cluster==i) %>%
      #   mutate(prof = c(50,150,250,350,450,550,650,750))
    }
    
   
    # incert$slope_min <- incert$slope_min
    # incert$slope_max <- incert$slope_max
    
   
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


    regression_att <- nls(formula = flux ~ a * (prof)^(-b), data = slopes, start = list(a=60, b= 0.8))
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
                              scale_y_continuous(name = expression(paste(Seasonal~accumulation,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,90)) +
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
                                scale_y_continuous(name = expression(paste(Accumulation~saisonnière,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,90)) +
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
    #############################################################################
    ############################ TRANSFERT EFFICIENCY ###########################
    #############################################################################
    #############################################################################
    
    transfert <- slopes %>%
      dplyr::mutate(transfert_eff = flux/flux[1])%>%
      dplyr::select(prof, transfert_eff)
    
    if(i!=5){
      print(paste("transfert à 1000m :", transfert$transfert_eff[8] * 100, "%"))}
    
    if(i==5){
      print(paste("transfert à 1000m :", transfert$transfert_eff[6] * 100, "%"))}


    #############################################################################
    ############################ IPOC PLOT ######################################
    #############################################################################
    #############################################################################
    
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
    
    scaling_factor_3 <- sec_axis_adjustement_factors(-c(0,600), c(1,30))
    
    I_plots[[i]] <- I_plots[[i]] +
      geom_line(aes(y = productive_layer, color = "y1"), size=0.5) + 
      scale_x_continuous(name = "relative month", expand = c(0.005,0)) +
      scale_y_continuous(expression(paste(iPOC,"  [",g~C~m^{-2},"]")), breaks =  c(1,5,10,15,20,25,30), sec.axis = sec_axis(trans = ~ {. - scaling_factor_3$adjust} / scaling_factor_3$diff *-1, name = "MLD [m]")) + 
      geom_line(aes(y = -mld * scaling_factor_3$diff + scaling_factor_3$adjust, color = "y2"), size = 0.5) + 
      scale_color_manual(values = c("y1" = "#CC0000", "y2" = "#3399FF")) +
      coord_cartesian(xlim = c(0,12), ylim = c(1,30)) + 
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
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_0_800_x)], 
                   y = min_0_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_0_800_x)], 
                   yend = max_0_800_y, 
                   color = "#66CC00", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_100_800_x)], 
                   y = min_100_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_100_800_x)], 
                   yend = max_100_800_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_200_800_x)], 
                   y = min_200_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_200_800_x)], 
                   yend = max_200_800_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_300_800_x)], 
                   y = min_300_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_300_800_x)], 
                   yend = max_300_800_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_400_800_x)], 
                   y = min_400_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_400_800_x)], 
                   yend = max_400_800_y, 
                   color = "#339900", size = 0.3) +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_500_800_x)], 
                   y = min_500_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_500_800_x)], 
                   yend = max_500_800_y, 
                   color = "#339900", size = 0.3)
    
    if(i!=5){
      I_plots[[i]] <- I_plots[[i]] +
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_600_800_x)], 
                   y = min_600_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_600_800_x)], 
                   yend = max_600_800_y, 
                   color = "#339900", size = 0.3)+
      geom_segment(x = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == min_700_800_x)], 
                   y = min_700_800_y, 
                   xend = I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == max_700_800_x)], 
                   yend = max_700_800_y, 
                   color = "#339900", size = 0.3)
    }
    
   I_values_per_layers_ready_to_plot$export <-  I_values_per_layers_ready_to_plot$depth_0_800 / I_values_per_layers_ready_to_plot$productive_layer
   print(plot(I_values_per_layers_ready_to_plot$dates,I_values_per_layers_ready_to_plot$export))
   
    # ###################
    # #############################################################################
    # ################## INSTANTANEOUS AND NORMALIZED FLUX ########################
    # #############################################################################
    # #############################################################################
    # 
    # climax <- climax_x
    # before_apex <- as.integer((apex_x + climax_x)/2)
    # before_apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == before_apex)])
    # after_apex <- apex_x + (apex_x - climax_x)
    # after_apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == after_apex)])
    # around_apex <- apex_x
    # around_apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == around_apex)])
    # 
    # instant_flux <- data.frame(depths = depth_nb, climax_Ez = NA, climax_flux = NA, climax_sd = NA,
    #                            before_Ez = NA, before_flux = NA, before_sd = NA,
    #                            around_Ez = NA, around_flux = NA, around_sd = NA,
    #                            after_Ez = NA, after_flux = NA, after_sd = NA)
    # couche <- filter(values_per_layers, layer_name == grep(x = depth_layers, pattern = "depth", value = TRUE)[1])
    # Ez <- couche$Ez_values_smth*10^3
    # 
    # instant_flux$climax_Ez[1] <- Ez[climax]
    # instant_flux$before_Ez[1] <- Ez[before_apex]
    # instant_flux$around_Ez[1] <- Ez[around_apex]
    # instant_flux$after_Ez[1] <- Ez[after_apex]
    # 
    # instant_flux$climax_flux[1] <- 1
    # instant_flux$before_flux[1] <- 1
    # instant_flux$around_flux[1] <- 1
    # instant_flux$after_flux[1] <- 1
    # 
    # instant_flux$climax_sd[1] = sd(Ez[climax])
    # instant_flux$before_sd[1] = sd(Ez[before_apex])
    # instant_flux$around_sd[1] = sd(Ez[around_apex])
    # instant_flux$after_sd[1] = sd(Ez[after_apex])
    # 
    # 
    # for (j in 2:6) {
    #   
    #   couche <- filter(values_per_layers, layer_name == grep(x = depth_layers, pattern = "depth", value = TRUE)[j])
    #   smth_I <- couche$smoothed_I
    #   Ez <- couche$Ez_values_smth*10^3
    #   
    #   instant_flux$climax_Ez[j] <- Ez[climax]
    #   instant_flux$climax_flux[j] <- instant_flux$climax_Ez[j]/instant_flux$climax_Ez[1]
    #   instant_flux$climax_sd[j] = sd(Ez[climax])
    #   
    #   instant_flux$before_Ez[j] <- Ez[before_apex]
    #   instant_flux$before_flux[j] <- instant_flux$before_Ez[j]/instant_flux$before_Ez[1]
    #   instant_flux$before_sd[j] = sd(Ez[before_apex])
    #   
    #   instant_flux$around_Ez[j] <- Ez[around_apex]
    #   instant_flux$around_flux[j] <- instant_flux$around_Ez[j]/instant_flux$around_Ez[1]
    #   instant_flux$around_sd[j] = sd(Ez[around_apex])
    #   
    #   instant_flux$after_Ez[j] <- Ez[after_apex]
    #   instant_flux$after_flux[j] <-  instant_flux$after_Ez[j]/instant_flux$after_Ez[1]
    #   instant_flux$after_sd[j] = sd(Ez[after_apex])
    #   
    # }
    # 
    # # 
    # # Instant_norm_flux_plots[[i]] <- ggplot(data = instant_flux, aes(x = depths)) +
    # #   
    # #   scale_x_reverse(name = "Depth below Zp [m]", expand = c(0,0), limits = c(800,-5)) +
    # #   scale_y_continuous(name = expression(paste(Normalized~iPOC~flux~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0,0),limits = c(-3,2)) +
    # #   coord_flip() +
    # #   
    # #   geom_line(aes(y = climax_flux), color = "#FF9900", size = 0.7, alpha = 0.8) +
    # #   geom_point(aes(y = climax_flux), shape=21, fill = "#FF9900", size=1.2, color="black") +
    # #   geom_errorbar( aes(ymin = climax_flux-climax_sd, ymax = climax_flux+climax_sd),width = 0.2) +
    # #   
    # #   geom_line(aes(y = before_flux), color = "#99CCCC", size = 0.7, alpha = 0.8) +
    # #   geom_point(aes(y = before_flux), shape=21, fill = "#99CCCC", size=1.2, color="black") +
    # #   geom_errorbar( aes(ymin = before_flux-before_sd, ymax = before_flux+before_sd),width = 0.2) +
    # #   
    # #   geom_line(aes(y = around_flux), color = "#FF3300", size = 0.7, alpha = 0.8) +
    # #   geom_point(aes(y = around_flux), shape=21, fill = "#FF3300", size=1.2, color="black") +
    # #   geom_errorbar( aes(ymin = around_flux-around_sd, ymax = around_flux+around_sd),width = 0.2) +
    # #   
    # #   geom_line(aes(y = after_flux), color = "grey", size = 0.7, alpha = 0.8) +
    # #   geom_point(aes(y = after_flux), shape=21, fill = "grey", size=1.2, color="black") +
    # #   geom_errorbar( aes(ymin = after_flux-after_sd, ymax = after_flux+after_sd),width = 0.2) +
    # #   
    # #   geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.4) +
    # #   
    # #   theme_bw() +
    # #   theme(plot.title = element_text(size=10, face="bold.italic", hjust = 0.5)) +
    # #   ggtitle(paste("cluster",i,sep=" ")) 
    # 
    # instant_flux <- instant_flux[-6,]
    # Instant_flux_plots[[i]] <- ggplot(data = instant_flux, aes(x = depths)) +
    #   
    #   scale_x_reverse(name = "Depth below Zp [m]", expand = c(0,0), limits = c(700,-5)) +
    #   scale_y_continuous(name = expression(paste(Instantaneous~iPOC~flux~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0,0),limits = c(-35,115)) +
    #   coord_flip() +
    #   
    #   geom_line(aes(y = climax_Ez), color = "#FF9900", size = 0.7, alpha = 0.8) +
    #   geom_point(aes(y = climax_Ez), shape=21, fill = "#FF9900", size=1.2, color="black") +
    #   geom_errorbar( aes(ymin = climax_Ez-climax_sd, ymax = climax_Ez+climax_sd),width = 0.2,color = "#FF9900") +
    #   
    #   geom_line(aes(y = before_Ez), color = "#99CCCC", size = 0.7, alpha = 0.8) +
    #   geom_point(aes(y = before_Ez), shape=21, fill = "#99CCCC", size=1.2, color="black") +
    #   geom_errorbar( aes(ymin = before_Ez-before_sd, ymax = before_Ez+before_sd),width = 0.2,color = "#99CCCC") +
    #   
    #   geom_line(aes(y = around_Ez), color = "#FF3300", size = 0.7, alpha = 0.8) +
    #   geom_point(aes(y = around_Ez), shape=21, fill = "#FF3300", size=1.2, color="black") +
    #   geom_errorbar( aes(ymin = around_Ez-around_sd, ymax = around_Ez+around_sd),width = 0.2, color = "#FF3300") +
    #   
    #   geom_line(aes(y = after_Ez), color = "grey", size = 0.7, alpha = 0.8) +
    #   geom_point(aes(y = after_Ez), shape=21, fill = "grey", size=1.2, color="black") +
    #   geom_errorbar( aes(ymin = after_Ez-after_sd, ymax = after_Ez+after_sd),width = 0.2, color = "grey") +
    #   
    #   geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.4) +
    #   
    #   theme_bw() +
    #   theme(plot.title = element_text(size=10, face="bold.italic", hjust = 0.5),
    #         legend.position = "left") +
    #   ggtitle(paste("cluster",i,sep=" "))
    # 
    
    # #############################################################################
    # ############################ EZ PLOT ########################################
    # #############################################################################
    # #############################################################################
    # 
    # Ez_values_per_layers_ready_to_plot <- values_per_layers %>% 
    #   pivot_wider(id_cols = c(month, zp), names_from = "layer_name",values_from = "Ez_values_smth") %>% 
    #   .[-1,]
    # 
    # Ez_plots[[i]] <- ggplot(data = Ez_values_per_layers_ready_to_plot*1000, aes(x = month/1000)) #Ez in mg POC/m2
    # # 
    # for (j in 1:length(depth_nb) ) {
    #   Ez_plots[[i]] <- Ez_plots[[i]] +
    #     geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = colors[j], size = 0.5)
    #   #geom_point(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = colors[j], size = 1.5)
    # }
    # 
    # Ez_plots[[i]] <- Ez_plots[[i]] +
    #   
    #   geom_line(aes(y = productive_layer), color = "red", size = 0.7, alpha = 0.8) +
    #   #geom_line(aes(y = whole_column), color = "#FFCC33", size = 1) +
    #   #geom_point(aes(y = productive_layer), color = "red", size = 1.5)+
    #   scale_x_continuous(name = "relative month", expand = c(0,0)) +
    #   scale_y_continuous(name = expression(paste(E[Z],"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0)) +
    #   coord_cartesian(xlim = c(0,12), ylim = c(-150,200)) + 
    #   
    #   theme_bw() +
    #   
    #   geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
    #   geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#FFCC33", size = 0.8) +
    #   geom_vline(xintercept = climax_x_month, linetype = "dotted", color = "#FF9900", size = 0.8) +
    #   geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF3300", size = 0.8) +
    #   geom_vline(xintercept = after_apex_x_month, linetype = "dotted", color = "grey", size = 0.8) +
    #   geom_vline(xintercept = before_apex_x_month, linetype = "dotted", color = "#99CCCC", size = 0.8) +
    #   
    #   theme (text = element_text(size=25, colour = "black"),
    #          axis.text.x = element_text(size = 8, colour = "black"), 
    #          axis.text.y = element_text(size = 8, colour = "black"),
    #          axis.title.x = element_text(size = 8),
    #          axis.title.y = element_text(size = 8),
    #          plot.title = element_text(size=10, face="bold.italic", hjust = 0.5)) +
    #   ggtitle(paste("cluster",i,sep=" "))
    
    #############################################################################
    ########################## SINKING SPEED ####################################
    #############################################################################
    #############################################################################
    #
    # if(i == 1 | i == 5){
    #   
    #   maximums_y <- summarise(I_values_per_layers_ready_to_plot, depth_0_150=max(depth_0_150),
    #                           depth_150_300=max(depth_150_300), depth_300_450=max(depth_300_450),
    #                           depth_450_600=max(depth_450_600), depth_600_750=max(depth_600_750))%>%
    #     t()
    #   
    #   max_1 <- match(max(I_values_per_layers_ready_to_plot$depth_0_150, na.rm = T), I_values_per_layers_ready_to_plot$depth_0_150)
    #   max_2 <- match(max(I_values_per_layers_ready_to_plot$depth_150_300, na.rm = T), I_values_per_layers_ready_to_plot$depth_150_300)
    #   max_3 <- match(max(I_values_per_layers_ready_to_plot$depth_300_450, na.rm = T), I_values_per_layers_ready_to_plot$depth_300_450)
    #   max_4 <- match(max(I_values_per_layers_ready_to_plot$depth_450_600, na.rm = T), I_values_per_layers_ready_to_plot$depth_450_600)
    #   max_5 <- match(max(I_values_per_layers_ready_to_plot$depth_600_750, na.rm = T), I_values_per_layers_ready_to_plot$depth_600_750)
    #   
    #   max_1_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_1)])
    #   max_2_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_2)])
    #   max_3_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_3)])
    #   max_4_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_4)])
    #   max_5_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_5)])
    #   
    #   dates_max <- c(max_1, max_2, max_3, max_4, max_5)
    #   dates_max_month <- c(max_1_month, max_2_month, max_3_month, max_4_month, max_5_month)
    #   
    #   maximums_300 <-  data.frame(depths = as.numeric(depth_nb[1:3])) %>%
    #     cbind(dates_max = dates_max[1:3])
    #   
    #   maximums_all <-  data.frame(depths_all = as.numeric(depth_nb[-6])) %>%
    #     cbind(dates_max_all = dates_max)
    #   
    #   sink_speed_regression <- lm(data = maximums_300, formula = (dates_max-dates_max[1]) ~ depths +0)
    #   
    #   sinking_speed_plots[[i]] <- ggplot(data = maximums_all, aes(x = dates_max_all)) +
    #     geom_point(aes(y = dates_max_all, x = depths_all), size = 2, color = colors[-6]) +
    #     geom_abline(intercept = as.numeric(dates_max[1]), slope = as.numeric(-sink_speed_regression$coefficients[1]), color = "red") +
    #     scale_y_continuous(name = "relative Julian Day", expand = c(0,0), limits = c(0,365)) +
    #     scale_x_reverse(name = "layers below Zp [m] (y:1000m)", expand = c(0,0), limits = c(800,-10)) +
    #     coord_flip() +
    #     annotation_custom(grobTree(textGrob(paste("speed =",format(as.numeric(1/sink_speed_regression$coefficients[1]), scientific = F, digits = 2)), x=0.68,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) +
    #     annotation_custom(grobTree(textGrob(expression(paste("  [",m~d^{-1},"]")), x=0.84,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) +
    #     theme_bw()
    #   
    # }
    
    
    ############################################################################
    ########################### RSUB PLOT ######################################
    ############################################################################
    ############################################################################
    
    DOXY_profile_data <- plotData_DOXY %>%
      filter(PARAM == "AOU", cluster == i) %>%
      left_join(zp_data %>% dplyr::select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(mld_data %>% dplyr::select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      left_join(plotData_SIG %>% dplyr::select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
      mutate(depths_below_zp = NA)
    DOXY_profile_data$zp <- round(DOXY_profile_data$zp)
    
    
    depths = DOXY_profile_data$PRES
    dates = DOXY_profile_data$jDay
    month = DOXY_profile_data$month
    zp = DOXY_profile_data$zp
    mld = DOXY_profile_data$mld
    depth_bin = 100
    SIG <- DOXY_profile_data$SIG +1000
    values <- DOXY_profile_data$value * SIG
    
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
        # (mean(df$values[index]))
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
  
  depth_0_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "0")
  iO2_0_800 <- depth_0_800$smoothed_I
  
  depth_100_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "100")
  iO2_100_800 <- depth_100_800$smoothed_I
  
  depth_200_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "200")
  iO2_200_800 <- depth_200_800$smoothed_I
  
  depth_300_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "300")
  iO2_300_800 <- depth_300_800$smoothed_I
  
  depth_400_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "400")
  iO2_400_800 <- depth_400_800$smoothed_I
  
  depth_500_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "500")
  iO2_500_800 <- depth_500_800$smoothed_I
  
  depth_600_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "600")
  iO2_600_800 <- depth_600_800$smoothed_I
  
  depth_700_800 <- values_per_layers_O2 %>%
    filter(DEPTH == "700")
  iO2_700_800 <- depth_700_800$smoothed_I
  
  
  ######### 0
  min_0_800_x <- match(min(iO2_0_800),iO2_0_800)
  max_0_800_x <- match(max(iO2_0_800),iO2_0_800)
  min_0_800_y <- iO2_0_800[min_0_800_x]
  max_0_800_y <- iO2_0_800[max_0_800_x]
  if(min_0_800_x<max_0_800_x){
    slope_O2_0_800 <- (min_0_800_y - max_0_800_y) / (min_0_800_x - max_0_800_x)
  }
  if(min_0_800_x>max_0_800_x){
    slope_O2_0_800 <- (max_0_800_y - min_0_800_y) / (365-min_0_800_x + max_0_800_x)
  }
  ######### 150
  min_100_800_x <- match(min(iO2_100_800),iO2_100_800)
  min_100_800_y <- iO2_100_800[min_100_800_x]
  max_100_800_x <- match(max(iO2_100_800),iO2_100_800)
  max_100_800_y <- iO2_100_800[max_100_800_x]
  if(min_100_800_x<max_100_800_x){
    slope_O2_100_800 <- (min_100_800_y - max_100_800_y) / (min_100_800_x - max_100_800_x)
  }
  if(min_100_800_x>max_100_800_x){
    slope_O2_100_800 <- (max_100_800_y - min_100_800_y) / (365-min_100_800_x + max_100_800_x)
  }
  ######### 300
  min_200_800_x <- match(min(iO2_200_800),iO2_200_800)
  min_200_800_y <- iO2_200_800[min_200_800_x]
  max_200_800_x <- match(max(iO2_200_800),iO2_200_800)
  max_200_800_y <- iO2_200_800[max_200_800_x]
  if(min_200_800_x<max_200_800_x){
    slope_O2_200_800 <- (min_200_800_y - max_200_800_y) / (min_200_800_x - max_200_800_x)
  }
  if(min_200_800_x>max_200_800_x){
    slope_O2_200_800 <- (max_200_800_y - min_200_800_y) / (365-min_200_800_x + max_200_800_x)
  }
  ######### 450
  min_300_800_x <- match(min(iO2_300_800),iO2_300_800)
  min_300_800_y <- iO2_300_800[min_300_800_x]
  max_300_800_x <- match(max(iO2_300_800),iO2_300_800)
  max_300_800_y <- iO2_300_800[max_300_800_x]
  if(min_300_800_x<max_300_800_x){
    slope_O2_300_800 <- (min_300_800_y - max_300_800_y) / (min_300_800_x - max_300_800_x)
  }
  if(min_300_800_x>max_300_800_x){
    slope_O2_300_800 <- (max_300_800_y - min_300_800_y) / (365-min_300_800_x + max_300_800_x)
  }
  ######### 600
  min_400_800_x <- match(min(iO2_400_800),iO2_400_800)
  min_400_800_y <- iO2_400_800[min_400_800_x]
  max_400_800_x <- match(max(iO2_400_800),iO2_400_800)
  max_400_800_y <- iO2_400_800[max_400_800_x]
  if(min_400_800_x<max_400_800_x){
    slope_O2_400_800 <- (min_400_800_y - max_400_800_y) / (min_400_800_x - max_400_800_x)
  }
  if(min_400_800_x>max_400_800_x){
    slope_O2_400_800 <- (max_400_800_y - min_400_800_y) / (365-min_400_800_x + max_400_800_x)
  }
  ######### 750
  min_500_800_x <- match(min(iO2_500_800),iO2_500_800)
  min_500_800_y <- iO2_500_800[min_500_800_x]
  max_500_800_x <- match(max(iO2_500_800),iO2_500_800)
  max_500_800_y <- iO2_500_800[max_500_800_x]
  if(min_500_800_x<max_500_800_x){
    slope_O2_500_800 <- (min_500_800_y - max_500_800_y) / (min_500_800_x - max_500_800_x)
  }
  if(min_500_800_x>max_500_800_x){
    slope_O2_500_800 <- (max_500_800_y - min_500_800_y) / (365-min_500_800_x + max_500_800_x)
  }
  
  
  if(i!=5){
      min_600_800_x <- match(min(iO2_600_800),iO2_600_800)
      min_600_800_y <- iO2_600_800[min_600_800_x]
      max_600_800_x <- match(max(iO2_600_800),iO2_600_800)
      max_600_800_y <- iO2_600_800[max_600_800_x]
      if(min_600_800_x<max_600_800_x){
        slope_O2_600_800 <- (min_600_800_y - max_600_800_y) / (min_600_800_x - max_600_800_x)
      }
      if(min_600_800_x>max_600_800_x){
        slope_O2_600_800 <- (max_600_800_y - min_600_800_y) / (365-min_600_800_x + max_600_800_x)
      }
      
      min_700_800_x <- match(min(iO2_700_800),iO2_700_800)
      min_700_800_y <- iO2_700_800[min_700_800_x]
      max_700_800_x <- match(max(iO2_700_800),iO2_700_800)
      max_700_800_y <- iO2_700_800[max_700_800_x]
      if(min_700_800_x<max_700_800_x){
        slope_O2_700_800 <- (min_700_800_y - max_700_800_y) / (min_700_800_x - max_700_800_x)
      }
      if(min_700_800_x>max_700_800_x){
        slope_O2_700_800 <- (max_700_800_y - min_700_800_y) / (365-min_700_800_x + max_700_800_x)
      }
      
      slopes_O2 <- data.frame(slope = c(slope_O2_0_800,slope_O2_100_800,slope_O2_200_800,
                                        slope_O2_300_800,slope_O2_400_800,slope_O2_500_800, 
                                        slope_O2_600_800,slope_O2_700_800))
        slopes_O2 <- slopes_O2 %>%mutate(i_depths = depth_nb)  
  }
  
  if(i==5){
    slopes_O2 <- data.frame(slope = c(slope_O2_0_800,slope_O2_100_800,slope_O2_200_800,
                                      slope_O2_300_800,slope_O2_400_800,slope_O2_500_800))%>%
      mutate(i_depths = depth_nb[-(7:8)])
  }
  
  slopes_O2$slope[1] <- (1/1.45)*slopes_O2$slope[1] *10^-3 *12     # x le quotient respiratoire pour passer de O2 en C, 
  slopes_O2$slope[2] <- (1/1.45)*slopes_O2$slope[2] *10^-3 *12  # x0.001 pour passer de µmol en mmol,
  slopes_O2$slope[3] <- (1/1.45)*slopes_O2$slope[3] *10^-3 *12    # x12 pour passer mmol de C en mg de C
  slopes_O2$slope[4] <- (1/1.45)*slopes_O2$slope[4] *10^-3 *12 
  slopes_O2$slope[5] <- (1/1.45)*slopes_O2$slope[5] *10^-3 *12 
  slopes_O2$slope[6] <- (1/1.45)*slopes_O2$slope[6] *10^-3 *12 
  
  if(i!=5){
    slopes_O2$slope[7] <- (1/1.45)*slopes_O2$slope[7] *10^-3 *12 
    slopes_O2$slope[8] <- (1/1.45)*slopes_O2$slope[8] *10^-3 *12
  }
  
 #####################################################################
  
  incertitudes <- data.frame(cluster=i, layer = depth_layers[3:10], slope_min = rep(NA,8), slope_max = rep(NA,8))
  layer_name <- paste(depth_seq, 800, sep = "_")
  
  for (l in 3:10){
    
    print(paste("couche", layer_name[l]))
    
    min_y <- get(paste("min_",layer_name[l],"_y", sep=""))
    max_y <- get(paste("max_",layer_name[l],"_y", sep=""))
    min_x <- get(paste("min_",layer_name[l],"_x", sep=""))
    max_x <- get(paste("max_",layer_name[l],"_x", sep=""))
    sd_tot <- incertitudes_tot_AOU(AOU_sd, window=30, i)[[l-2]]
    
    if(min_x<max_x){
      # pente minimum
      new_min_y <- min_y + sd_tot[min_x]
      new_max_y <- max_y - sd_tot[max_x]
      slope_min <- (new_max_y - new_min_y) / (max_x - min_x) 
      
      # pente maximum
      new_min_y <- min_y - sd_tot[min_x]
      new_max_y <- max_y + sd_tot[max_x]
      slope_max <- (new_max_y - new_min_y) / (max_x - min_x) 
    }
    
    if(min_x>max_x){
      # pente minimum
      new_min_y <- min_y + sd_tot[min_x]
      new_max_y <- max_y - sd_tot[max_x]
      slope_min <- (new_max_y - new_min_y) / ((365+max_x) - min_x) 
      
      # pente maximum
      new_min_y <- min_y - sd_tot[min_x]
      new_max_y <- max_y + sd_tot[max_x]
      slope_max <- (new_max_y - new_min_y) / ((365+max_x) - min_x)
    }
    
    incertitudes$slope_min[l-2] <- slope_min * (1/1.45) *10^-3 *12
    incertitudes$slope_max[l-2] <- slope_max * (1/1.45) *10^-3 *12
  }
  
  
  ########################################################################
  
  slopes_O2 <- slopes_O2 %>% mutate(slope_min = incertitudes$slope_min,
                                    slope_max = incertitudes$slope_max,
                                    diff = slope_max - slope)
  
  if(i!=5){
  slopes_O2 <- slopes_O2 %>%
    mutate(prof = c(50,150,250,350,450,550,650,750))
  # incert_O <- filter(incertitudes_O, cluster==i) %>%
  #   mutate(prof = c(50,150,250,350,450,550,650,750))
  }
  
  if(i==5){
    slopes_O2 <- slopes_O2 %>%
      mutate(prof = c(50,150,250,350,450,550))
    # incert_O <- filter(incertitudes_O, cluster==i) %>%
    #   mutate(prof = c(50,150,250,350,450,550))
  }
  
  # incert_O$slope_min <- incert_O$slope_min
  # incert_O$slope_max <- incert_O$slope_max
  
 
  
  if(i!=5){
  Rsub_plots[[i]] <- ggplot(data = slopes_O2, aes(x = prof, y = slope)) +
    geom_col(col="black", fill = colors) +
    geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
    # geom_point(shape=21, fill = colors, size=2) +
    # geom_errorbar(aes(ymin = Rsub-sd_tot, ymax = Rsub+sd_tot), width = 15) +
    scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(800,-10)) +
    scale_y_continuous(name = expression(paste(Respiration~de~C~intégrée~nette~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0.01,0.01), limits = c(-30,1100)) +
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
      scale_y_continuous(name = expression(paste(Respiration~de~C~intégrée~nette~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0.01,0.01), limits = c(-30,1100)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
      coord_flip() +
      theme_bw() +
      theme(plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
      ggtitle(names[i])
  }
    #############################################################################
    ############################ COMPARAISON ######################################
    #############################################################################
    #############################################################################        
    

    if(i==1){
      Rsub_Att <- bind_cols(att=slopes$flux, slope_min_C = slopes$slope_min, slope_max_C = slopes$slope_max, 
                            Rsub=slopes_O2$slope, slope_min_O = slopes_O2$slope_min, slope_max_O = slopes_O2$slope_max)
      regression_Rsub_Att[[i]] <- lm(data = Rsub_Att, formula = Rsub~att)

      # regression_Rsub_Att[[i]] <- nls(data = Rsub_Att, formula = log(Rsub)~ c + d *log(log(att)-h), start = list(c = 3, d = 1, h=0.9))
      # c[i] = summary(regression_Rsub_Att[[i]])$coefficients[1]
      # d[i] = summary(regression_Rsub_Att[[i]])$coefficients[2]
      # fonction[[i]] <- function(Rsub, c, d) c+d*log(log(Rsub))
    }

    if(i!=1 & i!=5){
      Rsub_Att_temp <- bind_cols(att=slopes$flux, slope_min_C = slopes$slope_min, slope_max_C = slopes$slope_max, 
                                 Rsub=slopes_O2$slope, slope_min_O = slopes_O2$slope_min, slope_max_O = slopes_O2$slope_max)
      regression_Rsub_Att[[i]] <- lm(data = Rsub_Att_temp, formula = Rsub~att)
      # regression_Rsub_Att[[i]] <- nls(data = Rsub_Att_temp, formula = log(Rsub)~ c + d *log(log(att)-h), start = list(c = 2, d = 1, h=0.2))
      # c[i] = summary(regression_Rsub_Att[[i]])$coefficients[1]
      # d[i] = summary(regression_Rsub_Att[[i]])$coefficients[2]
      # fonction[[i]] <- function(Rsub, c, d) c+d*log(log(Rsub))
      Rsub_Att <- bind_rows(Rsub_Att, Rsub_Att_temp)
    }

    if(i==5){
      Rsub_Att_temp <- bind_cols(att=c(slopes$flux,NA,NA), slope_min_C = c(slopes$slope_min,NA,NA), slope_max_C = c(slopes$slope_max,NA,NA),
                                 Rsub=c(slopes_O2$slope,NA,NA), slope_min_O = c(slopes_O2$slope_min,NA,NA), slope_max_O = c(slopes_O2$slope_max,NA,NA))
      regression_Rsub_Att[[i]] <- lm(data = Rsub_Att_temp, formula = Rsub~att)
      # regression_Rsub_Att[[i]] <- lm(data = Rsub_Att_temp, formula = log(Rsub)~log(log(att)))
      # c[i] = summary(regression_Rsub_Att[[i]])$coefficients[1]
      # d[i] = summary(regression_Rsub_Att[[i]])$coefficients[2]
      # fonction[[i]] <- function(Rsub, c, d) c+d*log(log(Rsub))
      Rsub_Att <- bind_rows(Rsub_Att, Rsub_Att_temp)
    }

    #  if(i==1){
    #    Rsub_Att <- bind_cols(att=transfert$transfert_eff, Rsub)
    #    regression_Rsub_Att[[i]] <- lm(data = Rsub_Att, formula = Rsub~att)
    #  }
    #
    #  if(i!=1){
    #    Rsub_Att_temp <- bind_cols(att=transfert$transfert_eff, Rsub)
    #    regression_Rsub_Att[[i]] <- lm(data = Rsub_Att_temp, formula = Rsub~att)
    #    Rsub_Att <- bind_rows(Rsub_Att, Rsub_Att_temp)
    #  }
    # print(paste("cluster",i,": y = ",slope = regression_Rsub_Att[[i]]$coefficients[[2]],"x +",regression_Rsub_Att[[i]]$coefficients[[1]]))
    
  }
  
  Rsub_Att <- Rsub_Att %>%
    mutate(cluster = c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8)),
           depths = rep(c("0","100","200","300","400","500","600","700"),5))
  
  Rsub_Att_plots <- ggplot(data = Rsub_Att, aes(x = att, y = Rsub, group=cluster)) +
    # geom_point(aes(fill=factor(cluster)),shape=21, size = 2) +
    geom_pointrange(aes(ymin = slope_min_O, ymax = slope_max_O, color=factor(cluster)), size = .2) +
    geom_pointrange(aes(xmin = slope_min_C, xmax = slope_max_C, color=factor(cluster)), size = .2) +
    
    scale_x_continuous(name = expression(paste(Accumulation~saisonnière~" [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0.05,0.05)) +
    #scale_x_continuous(name = expression(paste("Transfert efficiency [%]")), expand = c(0.05,0.05)) +
    scale_y_continuous(name = expression(paste(Respiration~nette~de~C~~"[",mg~C~m^{-2}~d^{-1},"]")), expand = c(0.01,0.01)) +
    
    
    geom_abline(intercept = regression_Rsub_Att[[1]]$coefficients[[1]], slope = regression_Rsub_Att[[1]]$coefficients[[2]], color = "#CC0000", linetype = "dashed") +
    annotation_custom(grobTree(textGrob(expression("(1)"~R^2~"="), x=0.75,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
    annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[1]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +

    geom_abline(intercept = regression_Rsub_Att[[2]]$coefficients[[1]], slope = regression_Rsub_Att[[2]]$coefficients[[2]], color = "#00CC00", linetype = "dashed",color = "#FF6666") +
    annotation_custom(grobTree(textGrob(expression("(2)"~R^2~"="), x=0.75,  y=0.86, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
    annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[2]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.86, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +

    geom_abline(intercept = regression_Rsub_Att[[3]]$coefficients[[1]], slope = regression_Rsub_Att[[3]]$coefficients[[2]], color = "#FFCC00", linetype = "dashed") +
    annotation_custom(grobTree(textGrob(expression("(3)"~R^2~"="), x=0.75,  y=0.82, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
    annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[3]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.82, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +

    geom_abline(intercept = regression_Rsub_Att[[4]]$coefficients[[1]], slope = regression_Rsub_Att[[4]]$coefficients[[2]], color = "#33CCFF", linetype = "dashed") +
    annotation_custom(grobTree(textGrob(expression("(4)"~R^2~"="), x=0.75,  y=0.78, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
    annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[4]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.78, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +

    geom_abline(intercept = regression_Rsub_Att[[5]]$coefficients[[1]], slope = regression_Rsub_Att[[5]]$coefficients[[2]], color ="#0066CC", linetype = "dashed") +
    annotation_custom(grobTree(textGrob(expression("(5)"~R^2~"="), x=0.75,  y=0.74, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
    annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[5]])$r.squared, scientific=F, digits = 3), x=0.89,  y=0.74, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +

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
    mutate(clust_name = c(rep("(1) Austral HCB",8),rep("(2) Austral MCB",8),rep("(3) STZ",8),rep("(4) PN",8),rep("(5) AN",8))) %>%
    filter(att>0) 
  
  Rsub_Att <- Rsub_Att%>%
    mutate(prof = c(rep(c("50","150","250","350","450","550","650","750"),4),c("50","150","250","350","450","550")))

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
##################
  # Ez_plots[[1]] <- Ez_plots[[1]] + theme(axis.title.x = element_blank())
  # Ez_plots[[2]] <- Ez_plots[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  # Ez_plots[[3]] <- Ez_plots[[3]] + theme(axis.title.y = element_blank())
  # #Ez_plots[[4]] <- Ez_plots[[4]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  # Ez_plots[[5]] <- Ez_plots[[5]] + theme(axis.title.y = element_blank())
###################
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
#####################" 
  Att_plots[[1]] <- Att_plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  Att_plots[[2]] <- Att_plots[[2]] +theme(axis.title.y.left = element_blank(), axis.title.x = element_blank(),
                                          axis.text.x = element_blank(), axis.text.y.left = element_blank())
  Att_plots[[3]] <- Att_plots[[3]] + theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank(), 
                                           axis.title.x = element_blank())
  Att_plots[[4]] <- Att_plots[[4]] + theme(axis.title.x = element_blank())
  Att_plots[[5]] <- Att_plots[[5]] + theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank())
  
  # Instant_flux_plots[[1]] <- Instant_flux_plots[[1]] + theme(axis.title.x = element_blank())
  # Instant_flux_plots[[2]] <- Instant_flux_plots[[2]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  # Instant_flux_plots[[3]] <- Instant_flux_plots[[3]] + theme(axis.title.y = element_blank())
  # #Instant_flux_plots[[4]] <- Instant_flux_plots[[4]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  # Instant_flux_plots[[5]] <- Instant_flux_plots[[5]] + theme(axis.title.y = element_blank())
  
  # Instant_norm_flux_plots[[1]] <- Instant_norm_flux_plots[[1]] + theme(axis.title.x = element_blank())
  # Instant_norm_flux_plots[[2]] <- Instant_norm_flux_plots[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  # Instant_norm_flux_plots[[3]] <- Instant_norm_flux_plots[[3]] + theme(axis.title.y = element_blank())
  # #Instant_norm_flux_plots[[4]] <- Instant_norm_flux_plots[[4]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
  # Instant_norm_flux_plots[[5]] <- Instant_norm_flux_plots[[5]] + theme(axis.title.y = element_blank())
  
  Rsub_plots[[1]] <- Rsub_plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  Rsub_plots[[2]] <- Rsub_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank(),
                                             axis.text.x = element_blank(), axis.text.y.left = element_blank())
  Rsub_plots[[3]] <- Rsub_plots[[3]] + theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank(), 
                                             axis.title.x = element_blank())
  Rsub_plots[[4]] <- Rsub_plots[[4]] + theme(axis.title.x = element_blank())
  Rsub_plots[[5]] <- Rsub_plots[[5]] + theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank())
  
  # Rsub_Att_plots[[1]] <- Rsub_Att_plots[[1]] + theme(axis.title.x = element_blank())
  # Rsub_Att_plots[[2]] <- Rsub_Att_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank())
  # Rsub_Att_plots[[3]] <- Rsub_Att_plots[[3]] + theme(axis.title.y.left = element_blank())
  # #Rsub_Att_plots[[4]] <- Rsub_Att_plots[[4]]
  # Rsub_plots[[5]] <- Rsub_Att_plots[[5]] + theme(axis.title.y.left = element_blank())
  
  finalPlot_I <- wrap_plots(I_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  # finalPlot_Ez <- wrap_plots(Ez_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  finalPlot_Rp <- wrap_plots(Rp_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  finalPlot_Att <- wrap_plots(Att_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  # finalPlot_Instant <- wrap_plots(Instant_flux_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  # finalPlot_norm_Instant <- wrap_plots(Instant_norm_flux_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  finalPlot_Rsub <- wrap_plots(Rsub_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  # finalPlot_Rsub_Att <- wrap_plots(Rsub_Att_plots,design=layout) #+ plot_annotation(tag_levels = "1")
  
  # finalPlot_Sink <- sinking_speed_plots[[1]] + sinking_speed_plots[[5]]
  #finalPlot_Sink <- wrap_plots(sinking_speed_plots,ncol = 2) #+ plot_annotation(tag_levels = "1")
  
  ggsave(plot = finalPlot_I,"figures/I_plots.png", width = 11.5, height = 5, dpi = 300)
  # ggsave(plot = finalPlot_Ez,"figures/Ez_plots.png", width = 11.5, height = 5, dpi = 300)
  ggsave(plot = finalPlot_Rp,"figures/Rp_plots.png", width = 11.5, height = 5, dpi = 300)
  ggsave(plot = finalPlot_Att,"figures/Att_plots.png", width = 11, height = 8, dpi = 300)
  # ggsave(plot = finalPlot_Instant,"figures/Instant_flux_plots.png", width = 11, height = 8, dpi = 300)
  # ggsave(plot = finalPlot_norm_Instant,"figures/Instant_flux_norm_plots.png", width = 11, height = 8, dpi = 300)
  ggsave(plot = finalPlot_Rsub,"figures/Rsub_plots.png", width = 11, height = 8, dpi = 300)
  # ggsave(plot = finalPlot_Sink,"figures/Sink_plots.png", width = 6.2, height = 4, dpi = 300)
  # ggsave(plot = finalPlot_Rsub_Att,"figures/Rsub_Att_plots.png", width = 11, height = 8, dpi = 300)
  ggsave(plot = Rsub_Att_plots,"figures/Rsub_Att_plots.png", width = 6, height = 4, dpi = 300)
  ggsave(plot = resp_plot,"figures/ratio_resp.png", width = 6, height = 4, dpi = 300)
} 
########################################################################




##############################################################################
# fonction Louis Terrats double axe ggplot
library(ggplot2)
library(grid)
library(patchwork)
library(gridExtra)


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



