
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
   Rsub_0_1000 <- vector()
   Rsub_150_1000 <- vector()
   Rsub_300_1000 <- vector()
   Rsub_450_1000 <- vector()
   Rsub_600_1000 <- vector()
   Rsub_750_1000 <- vector()
   sd_tot_0_1000 <- vector()
   sd_tot_150_1000 <- vector()
   sd_tot_300_1000 <- vector()
   sd_tot_450_1000 <- vector()
   sd_tot_600_1000 <- vector()
   sd_tot_750_1000 <- vector()
   regression_Rsub_Att <- list()
   max_mld <- vector()
   my_slope <- data.frame(depths_below_zp=rep(NA,954),term=rep(NA,954), estimate=rep(NA,954), std.error =rep(NA,954))
   
   for (i in 1:length(unique(profile_data$cluster))){
     
         print(paste("cluster", i))
         
         bbp_profile_data <- profile_data %>%
           filter(PARAM == "BBP700", cluster == i) %>%
           left_join(zp_data %>% select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
           left_join(mld_data %>% select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
           left_join(mask_data %>% select(jDay, PRES, POC_mask = value, cluster), by = c("jDay", "PRES", "cluster"))%>%
           left_join(POC_data %>% select(jDay, PRES, POC = value, cluster), by = c("jDay", "PRES", "cluster"))
         
         values <- bbp_profile_data$POC*10^-3
         
         depths = bbp_profile_data$PRES
         dates = bbp_profile_data$jDay
         month = bbp_profile_data$month
         zp = bbp_profile_data$zp
         mld = bbp_profile_data$mld
         POC_mask = bbp_profile_data$POC_mask
         depth_bin = 150
   
 #############################################################################
 ############################ PLOT MASK FOR METRIC ###########################
 #############################################################################
 #############################################################################
         values_iPOC_mask <- vector() %>% as.numeric()
         for (k in 1:length(unique(dates))){
           range <- which(dates==unique(dates)[k])
           values_iPOC_mask[k] <- sum(POC_mask[range]*10^-3)
         }
         smth_values_iPOC_mask <- smth(values_iPOC_mask, window = window, alpha = 2.5, method = "gaussian", tails = T)
         
         Ez_POC_mask <- c( NA, diff(smth_values_iPOC_mask) / {diff(unique(dates))})
         rp_values <- Ez_POC_mask/smth_values_iPOC_mask
         
         apex_x <- match(max(smth_values_iPOC_mask),smth_values_iPOC_mask) 
         apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == apex_x)])
         # onset_x <- match(min(smth_values_iPOC_mask[0:apex_x]),smth_values_iPOC_mask[0:apex_x])
         onset_x <- match(min(smth_values_iPOC_mask),smth_values_iPOC_mask)
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
         
         apex_y <- max(smth_values_iPOC_mask)
         climax_y <- smth_values_iPOC_mask[climax_x]
         onset_y <- min(smth_values_iPOC_mask[0:apex_x])
         
         metrics <- c(onset_x, climax_x, apex_x) %>% as.data.frame()
         write_csv(metrics, paste("data/metrics_cluster",i,".csv", sep=""))
         
         ###########
         #####
         scaling_factor_2 <- sec_axis_adjustement_factors(c(0,20), c(-0.02,0.03))
         
         ### accumulation rates
         
         rp_plot <- data.frame(rp = rp_values, month = unique(month), smth_iPOC_mask = smth_values_iPOC_mask) 
         
         Rp_plots[[i]] <- ggplot(data = rp_plot, aes(x = month)) +
           
           scale_x_continuous(name = "relative month", expand = c(0,0)) +
           
           scale_y_continuous(expression(paste(r," [",d^{-1},"]")), sec.axis = sec_axis(trans = ~ {. - scaling_factor_2$adjust} / scaling_factor_2$diff , name = expression(paste(iPOC,"  [",g~C~m^{-2},"]")))) + 
           
           geom_line(aes(y = rp), size = 0.5) +
           
           geom_line(aes(y = smth_iPOC_mask * scaling_factor_2$diff + scaling_factor_2$adjust), size = 0.8, color = "lightblue") + 
           
           coord_cartesian(xlim = c(0,12), ylim = c(-0.02,0.03)) + 
           
           geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.8) +
           
           geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#FFCC33", size = 0.8) +
           
           geom_vline(xintercept = climax_x_month, linetype = "dotted", color = "#FF9900", size = 0.8) +
           
           geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF3300", size = 0.8) +
           
           annotation_custom(grobTree(textGrob(expression(r[max]~"="), x=0.70,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
           
           annotation_custom(grobTree(textGrob(format(rp_plot$rp[climax_x], scientific = T, digits = 2), x=0.81,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
           
           annotation_custom(grobTree(textGrob(paste("onset = j", onset_x), x=0.70,  y=0.86, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
           annotation_custom(grobTree(textGrob(paste("climax = j", climax_x), x=0.70,  y=0.78, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
           annotation_custom(grobTree(textGrob(paste("apex = j", apex_x), x=0.70,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) + 
           
           theme_bw() +
           
           theme(text = element_text(size=25, colour = "black"),
                 axis.text.x = element_text(size = 8, colour = "black"), 
                 axis.text.y = element_text(size = 8, colour = "black"),
                 axis.title.x = element_text(size = 8),
                 axis.title.y = element_text(size = 8),
                 
                 axis.text.y.right = element_text(color = "black", size = 8), 
                 #axis.ticks.y.right = element_line(color = "black"),
                 #axis.line.y.right = element_line(color = "black"),
                 axis.title.y.right = element_text(color = "lightblue", size = 8),
                 
                 plot.title = element_text(size=10, face="bold.italic", hjust = 0.5)) +
           
           ggtitle(paste("cluster",i,sep=" "))
         
 #############################################################################
 ############################ ATTENUATION PLOT ###############################
 #############################################################################
 #############################################################################
         
         depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
         # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
         # 
         depth_layers <- ifelse(numbers_only(depth_seq), paste("depth", depth_seq, 800, sep = "_"), depth_seq)
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
         whole_column <- values_per_layers %>%
           filter(DEPTH == "whole_column")
      
         iPOC_whole_column <- whole_column$smoothed_I
         
         depth_0_1000 <- values_per_layers %>%
           filter(DEPTH == "0")
         iPOC_0_1000 <- depth_0_1000$smoothed_I
         
         depth_150_1000 <- values_per_layers %>%
           filter(DEPTH == "150")
         iPOC_150_1000 <- depth_150_1000$smoothed_I
         
         depth_300_1000 <- values_per_layers %>%
           filter(DEPTH == "300")
         iPOC_300_1000 <- depth_300_1000$smoothed_I
         
         depth_450_1000 <- values_per_layers %>%
           filter(DEPTH == "450")
         iPOC_450_1000 <- depth_450_1000$smoothed_I
         
         depth_600_1000 <- values_per_layers %>%
           filter(DEPTH == "600")
         iPOC_600_1000 <- depth_600_1000$smoothed_I
         
         depth_750_1000 <- values_per_layers %>%
           filter(DEPTH == "750")
         iPOC_750_1000 <- depth_750_1000$smoothed_I
         
         
         if(onset_x < apex_x){
             ######### 0
             #min_0_1000_x <- match(min(iPOC_0_1000[0:apex_x]),iPOC_0_1000[0:apex_x])
             min_0_1000_y <- iPOC_0_1000[onset_x]
             #max_0_1000_x <- match(max(iPOC_0_1000),iPOC_0_1000)
             max_0_1000_y <- iPOC_0_1000[apex_x]
             slope_0_1000 <- (max_0_1000_y - min_0_1000_y) / (apex_x - onset_x)
             ######### 150
             #min_150_1000_x <- match(min(iPOC_150_1000[0:apex_x]),iPOC_150_1000[0:apex_x])
             min_150_1000_y <- iPOC_150_1000[onset_x]
             #max_150_1000_x <- match(max(iPOC_150_1000),iPOC_150_1000)
             max_150_1000_y <- iPOC_150_1000[apex_x]
             slope_150_1000 <- (max_150_1000_y - min_150_1000_y) / (apex_x - onset_x)
             ######### 300
             #min_300_1000_x <- match(min(iPOC_300_1000[0:apex_x]),iPOC_300_1000[0:apex_x])
             min_300_1000_y <- iPOC_300_1000[onset_x]
             #max_300_1000_x <- match(max(iPOC_300_1000),iPOC_300_1000)
             max_300_1000_y <- iPOC_300_1000[apex_x]
             slope_300_1000 <- (max_300_1000_y - min_300_1000_y) / (apex_x - onset_x)
             ######### 450
             #min_450_1000_x <- match(min(iPOC_450_1000[0:apex_x]),iPOC_450_1000[0:apex_x])
             min_450_1000_y <- iPOC_450_1000[onset_x]
             #max_450_1000_x <- match(max(iPOC_450_1000),iPOC_450_1000)
             max_450_1000_y <- iPOC_450_1000[apex_x]
             slope_450_1000 <- (max_450_1000_y - min_450_1000_y) / (apex_x - onset_x)
             ######### 600
             #min_600_1000_x <- match(min(iPOC_600_1000[0:apex_x]),iPOC_600_1000[0:apex_x])
             min_600_1000_y <- iPOC_600_1000[onset_x]
             #max_600_1000_x <- match(max(iPOC_600_1000),iPOC_600_1000)
             max_600_1000_y <- iPOC_600_1000[apex_x]
             slope_600_1000 <- (max_600_1000_y - min_600_1000_y) / (apex_x - onset_x)
             ######### 750
             #min_750_1000_x <- match(min(iPOC_750_1000[0:apex_x]),iPOC_750_1000[0:apex_x])
             min_750_1000_y <- iPOC_750_1000[onset_x]
             #max_750_1000_x <- match(max(iPOC_750_1000),iPOC_750_1000)
             max_750_1000_y <- iPOC_750_1000[apex_x]
             slope_750_1000 <- (max_750_1000_y - min_750_1000_y) / (apex_x - onset_x)
         }
    
         if(onset_x > apex_x){
             ######### 0
             min_0_1000_x <- 0
             min_0_1000_y <- iPOC_0_1000[onset_x]
             max_0_1000_x <- 365-onset_x + apex_x
             max_0_1000_y <- iPOC_0_1000[apex_x]
             slope_0_1000 <- (max_0_1000_y - min_0_1000_y) / (max_0_1000_x)
             ######### 150
             min_150_1000_x <- 0
             min_150_1000_y <- iPOC_150_1000[onset_x]
             max_150_1000_x <- 365-onset_x + apex_x
             max_150_1000_y <- iPOC_150_1000[apex_x]
             slope_150_1000 <- (max_150_1000_y - min_150_1000_y) / (max_150_1000_x)
             ######### 300
             min_300_1000_x <- 0
             min_300_1000_y <- iPOC_300_1000[onset_x]
             max_300_1000_x <- 365-onset_x + apex_x
             max_300_1000_y <- iPOC_300_1000[apex_x]
             slope_300_1000 <- (max_300_1000_y - min_300_1000_y) / (max_300_1000_x)
             ######### 450
             min_450_1000_x <- 0
             min_450_1000_y <- iPOC_450_1000[onset_x]
             max_450_1000_x <- 365-onset_x + apex_x
             max_450_1000_y <- iPOC_450_1000[apex_x]
             slope_450_1000 <- (max_450_1000_y - min_450_1000_y) / (max_450_1000_x)
             ######### 600
             min_600_1000_x <- 0
             min_600_1000_y <- iPOC_600_1000[onset_x]
             max_600_1000_x <- 365-onset_x + apex_x
             max_600_1000_y <- iPOC_600_1000[apex_x]
             slope_600_1000 <- (max_600_1000_y - min_600_1000_y) / (max_600_1000_x)
             ######### 750
             min_750_1000_x <- 0
             min_750_1000_y <- iPOC_750_1000[onset_x]
             max_750_1000_x <- 365-onset_x + apex_x
             max_750_1000_y <- iPOC_750_1000[apex_x]
             slope_750_1000 <- (max_750_1000_y - min_750_1000_y) / (max_750_1000_x)
         }
              
         slopes <- data.frame(slope = c(slope_0_1000*10^3,slope_150_1000*10^3,slope_300_1000*10^3,slope_450_1000*10^3,slope_600_1000*10^3,slope_750_1000*10^3))%>%
           mutate(i_depths = depth_nb+1) 
         #################
      
         # model <- lm(log(slopes$slope)~log(slopes$i_depths))
         # b = model$coefficients[[2]]
         # Fz0 = model$coefficients[[1]]
         # # plot(log(slopes$slope)~log(slopes$i_depths))
         # # abline(model)
         # ###############
         # z = c(1:800)
         # regression_data <- data.frame(z, pred_slope = exp(Fz0)*z^b)
         
         regression_att <- nls(data = slopes, formula = slope ~ a * exp(-i_depths/b) + c, start = list(a=40, c=10, b= 150))
         a = summary(regression_att)$parameters[1]
         b = summary(regression_att)$parameters[3]
         c = summary(regression_att)$parameters[2]
         
         f <- function(i_depths, a, b, c) a * exp(-i_depths/b) + c
         
         Att_plots[[i]] <- print(ggplot(data = slopes, aes(x = i_depths, y = slope)) +
           geom_point(shape=21, fill = colors, size=1.2) +
           geom_function(fun = f, args = list(a=a, b=b, c=c)) +
           scale_x_reverse(name = "layers below Zp [m] (y:1000m)", expand = c(0,0), limits = c(800,0)) +
           scale_y_continuous(name = expression(paste(Seasonal~accumulation,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(-20,80)) +
           coord_flip() +
           theme_bw() +
           theme(axis.text = element_text(size = 10, colour = "black"),
                 axis.title = element_text(size = 10, colour = "black"),
                 axis.text.x=element_text(angle=0),
                 plot.title = element_text(size=10, face="bold.italic", hjust = 0.5)) +
           
           annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
           annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
           annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
           annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
           annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
           
           ggtitle(paste("cluster",i,sep=" ")))
   
 #############################################################################
 ############################ TRANSFERT EFFICIENCY ###########################
 #############################################################################
 #############################################################################
       
  transfert <- slopes %>%
           mutate(transfert_eff = slope/slope[1])
 
         
 #############################################################################
 ############################ IPOC PLOT ######################################
 #############################################################################
 #############################################################################
         
         I_values_per_layers_ready_to_plot <- values_per_layers %>% pivot_wider(id_cols = c(dates,month, zp, mld), 
                                                                                names_from = "layer_name",values_from = "smoothed_I")
         
         I_plots[[i]] <- ggplot(data = I_values_per_layers_ready_to_plot, aes(x = month)) +
           geom_area(aes_string(y = I_values_per_layers_ready_to_plot$whole_column), fill = "#FFCC33", color = "#FFCC33", alpha = 0.3)
         
         
         # geom_rect(data = flux_dates_df, aes(xmin = min_date, xmax = max_date, ymin = min_y, ymax = max_y), fill = "green3", alpha = 0.25)
         # geom_area(aes(y = whole_column), fill = "yellow", alpha = 0.5, color = "black", size = 0.15)
        
         
         for (j in seq_along(depth_nb) ) {
           I_plots[[i]] <- I_plots[[i]] + geom_area(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), fill = colors[j], color = colors[j]) +
             geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "black", size = 0.15) +
             geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "white", linetype = "dashed", size = 0.15)
         }
         
         scaling_factor_3 <- sec_axis_adjustement_factors(-c(0,600), c(0,40))
         
         I_plots[[i]] <- I_plots[[i]] +
           geom_line(aes(y = productive_layer, color = "y1"), size=0.5) + 
           scale_x_continuous(name = "relative month", expand = c(0,0)) +
           scale_y_continuous(expression(paste(iPOC,"  [",g~C~m^{-2},"]")), breaks =  c(0,5,10,15,20,25,30,35,40), sec.axis = sec_axis(trans = ~ {. - scaling_factor_3$adjust} / scaling_factor_3$diff *-1, name = "MLD (m)")) + 
           geom_line(aes(y = -mld * scaling_factor_3$diff + scaling_factor_3$adjust, color = "y2"), size = 0.5) + 
           scale_color_manual(values = c("y1" = "#CC0000", "y2" = "grey")) +
           coord_cartesian(xlim = c(0,12), ylim = c(0,40)) + 
           #geom_point(aes(x = apex_x, y = apex_y), shape=21, fill = "#FF3300", size=2, color="#FF3300") +
           geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF3300", size = 0.8) +
           #geom_point(aes(x = onset_x, y = onset_y), shape=21, fill = "#FFCC33", size=2, color="#FFCC33") +
           geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#FFCC33", size = 0.8) +
           #geom_point(aes(x = climax_x, y = climax_y), shape=21, fill = "#FF9900", size=2, color="#FF9900") +
           geom_vline(xintercept = climax_x_month, linetype = "dotted", color = "#FF9900", size = 0.8) +
           annotation_custom(grobTree(textGrob(expression(seasonal~E[0]~"="), x=0.63,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) + 
           annotation_custom(grobTree(textGrob(format(slope_0_1000, scientific = T, digits = 2), x=0.84,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) + 
           guides(colour = FALSE) +
           
           theme(text = element_text(size=8, colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_rect(linetype = "solid", fill = NA),
                 axis.text = element_text(size = 8, colour = "black"),
                 axis.title = element_text(size = 8, colour = "black"),
                 axis.text.x=element_text(angle=0),
                 
                 axis.text.y.left = element_text(color = "black"), 
                 axis.ticks.y.left = element_line(color = "black"),
                 axis.line.y.left = element_line(color = "black"),
                 axis.title.y.left = element_text(color = "black"),
                 
                 axis.text.y.right = element_text(color = "grey"), 
                 axis.ticks.y.right = element_line(color = "grey"),
                 axis.line.y.right = element_line(color = "grey"),
                 axis.title.y.right = element_text(color = "grey"),
                 plot.title = element_text(size=10, face="bold.italic", hjust = 0.5)) +
           ggtitle(paste("cluster",i,sep=" "))
         
         if(onset_x < climax_x){
           I_plots[[i]] <- I_plots[[i]] +
         geom_segment(x = onset_x_month, y = min_0_1000_y, xend = apex_x_month, yend = max_0_1000_y, color = "#66CCFF", size = 0.5)
         }
         
         if(onset_x > climax_x){
           y = min_0_1000_y + ((365-onset_x) * slope_0_1000)
           
           I_plots[[i]] <- I_plots[[i]] +
             geom_segment(x = onset_x_month, y = min_0_1000_y, xend = 12, yend = y, color = "#66CCFF", size = 0.5) +
             geom_segment(x = 0, y = y, xend = apex_x_month, yend = max_0_1000_y, color = "#66CCFF", size = 0.5)
         }
         ###################
 #############################################################################
 ################## INSTANTANEOUS AND NORMALIZED FLUX ########################
 #############################################################################
 #############################################################################

         climax <- climax_x
         before_apex <- as.integer((apex_x + climax_x)/2)
         before_apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == before_apex)])
         after_apex <- apex_x + (apex_x - climax_x)
         after_apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == after_apex)])
         around_apex <- apex_x
         around_apex_x_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == around_apex)])
      
         instant_flux <- data.frame(depths = depth_nb, climax_Ez = NA, climax_flux = NA, climax_sd = NA,
                                     before_Ez = NA, before_flux = NA, before_sd = NA,
                                     around_Ez = NA, around_flux = NA, around_sd = NA,
                                     after_Ez = NA, after_flux = NA, after_sd = NA)
         couche <- filter(values_per_layers, layer_name == grep(x = depth_layers, pattern = "depth", value = TRUE)[1])
         Ez <- couche$Ez_values_smth*10^3
         
         instant_flux$climax_Ez[1] <- Ez[climax]
         instant_flux$before_Ez[1] <- Ez[before_apex]
         instant_flux$around_Ez[1] <- Ez[around_apex]
         instant_flux$after_Ez[1] <- Ez[after_apex]
         
         instant_flux$climax_flux[1] <- 1
         instant_flux$before_flux[1] <- 1
         instant_flux$around_flux[1] <- 1
         instant_flux$after_flux[1] <- 1
         
         instant_flux$climax_sd[1] = sd(Ez[climax])
         instant_flux$before_sd[1] = sd(Ez[before_apex])
         instant_flux$around_sd[1] = sd(Ez[around_apex])
         instant_flux$after_sd[1] = sd(Ez[after_apex])
         
         
         for (j in 2:6) {
      
          couche <- filter(values_per_layers, layer_name == grep(x = depth_layers, pattern = "depth", value = TRUE)[j])
          smth_I <- couche$smoothed_I
          Ez <- couche$Ez_values_smth*10^3
       
          instant_flux$climax_Ez[j] <- Ez[climax]
          instant_flux$climax_flux[j] <- instant_flux$climax_Ez[j]/instant_flux$climax_Ez[1]
          instant_flux$climax_sd[j] = sd(Ez[climax])
             
          instant_flux$before_Ez[j] <- Ez[before_apex]
          instant_flux$before_flux[j] <- instant_flux$before_Ez[j]/instant_flux$before_Ez[1]
          instant_flux$before_sd[j] = sd(Ez[before_apex])
          
          instant_flux$around_Ez[j] <- Ez[around_apex]
          instant_flux$around_flux[j] <- instant_flux$around_Ez[j]/instant_flux$around_Ez[1]
          instant_flux$around_sd[j] = sd(Ez[around_apex])
          
          instant_flux$after_Ez[j] <- Ez[after_apex]
          instant_flux$after_flux[j] <-  instant_flux$after_Ez[j]/instant_flux$after_Ez[1]
          instant_flux$after_sd[j] = sd(Ez[after_apex])
          
         }
         
         
    Instant_norm_flux_plots[[i]] <- ggplot(data = instant_flux, aes(x = depths)) +
           
           scale_x_reverse(name = "Depth below Zp [m]", expand = c(0,0), limits = c(800,-5)) +
           scale_y_continuous(name = expression(paste(Normalized~iPOC~flux~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0,0),limits = c(-3,2)) +
           coord_flip() +
           
           geom_line(aes(y = climax_flux), color = "#FF9900", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = climax_flux), shape=21, fill = "#FF9900", size=1.2, color="black") +
           geom_errorbar( aes(ymin = climax_flux-climax_sd, ymax = climax_flux+climax_sd),width = 0.2) +
           
           geom_line(aes(y = before_flux), color = "#99CCCC", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = before_flux), shape=21, fill = "#99CCCC", size=1.2, color="black") +
           geom_errorbar( aes(ymin = before_flux-before_sd, ymax = before_flux+before_sd),width = 0.2) +
           
           geom_line(aes(y = around_flux), color = "#FF3300", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = around_flux), shape=21, fill = "#FF3300", size=1.2, color="black") +
           geom_errorbar( aes(ymin = around_flux-around_sd, ymax = around_flux+around_sd),width = 0.2) +
           
           geom_line(aes(y = after_flux), color = "grey", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = after_flux), shape=21, fill = "grey", size=1.2, color="black") +
           geom_errorbar( aes(ymin = after_flux-after_sd, ymax = after_flux+after_sd),width = 0.2) +
           
           geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.4) +
      
           theme_bw() +
           theme(plot.title = element_text(size=10, face="bold.italic", hjust = 0.5)) +
           ggtitle(paste("cluster",i,sep=" ")) 
         
         
    Instant_flux_plots[[i]] <- ggplot(data = instant_flux, aes(x = depths)) +
           
           scale_x_reverse(name = "Depth below Zp [m]", expand = c(0,0), limits = c(800,-5)) +
           scale_y_continuous(name = expression(paste(Instantaneous~iPOC~flux~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0,0),limits = c(-100,250)) +
           coord_flip() +
           
           geom_line(aes(y = climax_Ez), color = "#FF9900", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = climax_Ez), shape=21, fill = "#FF9900", size=1.2, color="black") +
           geom_errorbar( aes(ymin = climax_Ez-climax_sd, ymax = climax_Ez+climax_sd),width = 0.2,color = "#FF9900") +
           
           geom_line(aes(y = before_Ez), color = "#99CCCC", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = before_Ez), shape=21, fill = "#99CCCC", size=1.2, color="black") +
           geom_errorbar( aes(ymin = before_Ez-before_sd, ymax = before_Ez+before_sd),width = 0.2,color = "#99CCCC") +
           
           geom_line(aes(y = around_Ez), color = "#FF3300", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = around_Ez), shape=21, fill = "#FF3300", size=1.2, color="black") +
           geom_errorbar( aes(ymin = around_Ez-around_sd, ymax = around_Ez+around_sd),width = 0.2, color = "#FF3300") +
           
           geom_line(aes(y = after_Ez), color = "grey", size = 0.7, alpha = 0.8) +
           geom_point(aes(y = after_Ez), shape=21, fill = "grey", size=1.2, color="black") +
           geom_errorbar( aes(ymin = after_Ez-after_sd, ymax = after_Ez+after_sd),width = 0.2, color = "grey") +
           
           geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.4) +
           
           theme_bw() +
           theme(plot.title = element_text(size=10, face="bold.italic", hjust = 0.5),
                 legend.position = "left") +
           ggtitle(paste("cluster",i,sep=" "))
         
         
#############################################################################
############################ EZ PLOT ########################################
#############################################################################
#############################################################################
         
         Ez_values_per_layers_ready_to_plot <- values_per_layers %>% 
           pivot_wider(id_cols = c(month, zp), names_from = "layer_name",values_from = "Ez_values_smth") %>% 
           .[-1,]
         
         Ez_plots[[i]] <- ggplot(data = Ez_values_per_layers_ready_to_plot*1000, aes(x = month/1000)) #Ez in mg POC/m2
         # 
         for (j in 1:length(depth_nb) ) {
           Ez_plots[[i]] <- Ez_plots[[i]] +
             geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = colors[j], size = 0.5)
             #geom_point(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = colors[j], size = 1.5)
         }
      
         Ez_plots[[i]] <- Ez_plots[[i]] +
           
           geom_line(aes(y = productive_layer), color = "red", size = 0.7, alpha = 0.8) +
           #geom_line(aes(y = whole_column), color = "#FFCC33", size = 1) +
           #geom_point(aes(y = productive_layer), color = "red", size = 1.5)+
           scale_x_continuous(name = "relative month", expand = c(0,0)) +
           scale_y_continuous(name = expression(paste(E[Z],"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0)) +
           coord_cartesian(xlim = c(0,12), ylim = c(-250,350)) + 
      
           theme_bw() +
           
           geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
           geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#FFCC33", size = 0.8) +
           geom_vline(xintercept = climax_x_month, linetype = "dotted", color = "#FF9900", size = 0.8) +
           geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF3300", size = 0.8) +
           geom_vline(xintercept = after_apex_x_month, linetype = "dotted", color = "grey", size = 0.8) +
           geom_vline(xintercept = before_apex_x_month, linetype = "dotted", color = "#99CCCC", size = 0.8) +
           
           theme (text = element_text(size=25, colour = "black"),
                  axis.text.x = element_text(size = 8, colour = "black"), 
                  axis.text.y = element_text(size = 8, colour = "black"),
                  axis.title.x = element_text(size = 8),
                  axis.title.y = element_text(size = 8),
                  plot.title = element_text(size=10, face="bold.italic", hjust = 0.5)) +
           ggtitle(paste("cluster",i,sep=" "))
         
 #############################################################################
 ########################## SINKING SPEED ####################################
 #############################################################################
 #############################################################################
 #
         if(i == 1 | i == 5){
           
         maximums_y <- summarise(I_values_per_layers_ready_to_plot, depth_0_1000=max(depth_0_1000),
                               depth_150_1000=max(depth_150_1000), depth_300_1000=max(depth_300_1000),
                               depth_450_1000=max(depth_450_1000), depth_600_1000=max(depth_600_1000),
                               depth_750_1000=max(depth_750_1000))%>%
                               t()
          
         max_1 <- match(max(I_values_per_layers_ready_to_plot$depth_0_1000, na.rm = T), I_values_per_layers_ready_to_plot$depth_0_1000)
         max_2 <- match(max(I_values_per_layers_ready_to_plot$depth_150_1000, na.rm = T), I_values_per_layers_ready_to_plot$depth_150_1000)
         max_3 <- match(max(I_values_per_layers_ready_to_plot$depth_300_1000, na.rm = T), I_values_per_layers_ready_to_plot$depth_300_1000)
         max_4 <- match(max(I_values_per_layers_ready_to_plot$depth_450_1000, na.rm = T), I_values_per_layers_ready_to_plot$depth_450_1000)
         max_5 <- match(max(I_values_per_layers_ready_to_plot$depth_600_1000, na.rm = T), I_values_per_layers_ready_to_plot$depth_600_1000)
         max_6 <- match(max(I_values_per_layers_ready_to_plot$depth_750_1000, na.rm = T), I_values_per_layers_ready_to_plot$depth_750_1000)
         
         max_1_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_1)])
         max_2_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_2)])
         max_3_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_3)])
         max_4_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_4)])
         max_5_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_5)])
         max_6_month <- unique(bbp_profile_data$month[which(bbp_profile_data$jDay == max_6)])
         
         dates_max <- c(max_1, max_2, max_3, max_4, max_5, max_6)
         dates_max_month <- c(max_1_month, max_2_month, max_3_month, max_4_month, max_5_month, max_6_month)
         
         maximums_450 <-  data.frame(depths = as.numeric(depth_nb[1:4])) %>%
           cbind(dates_max = dates_max[1:4])
         
         maximums_all <-  data.frame(depths_all = as.numeric(depth_nb)) %>%
           cbind(dates_max_all = dates_max)

         sink_speed_regression <- lm(data = maximums_450, formula = (dates_max-dates_max[1]) ~ depths +0)

         sinking_speed_plots[[i]] <- ggplot(data = maximums_all, aes(x = dates_max)) +
           geom_point(aes(y = dates_max_all, x = depths_all), size = 2, color = colors) +
           geom_abline(intercept = as.numeric(dates_max[1]), slope = as.numeric(-sink_speed_regression$coefficients[1]), color = "red") +
           scale_y_continuous(name = "relative Julian Day", expand = c(0,0), limits = c(0,365)) +
           scale_x_reverse(name = "layers below Zp [m] (y:1000m)", expand = c(0,0), limits = c(800,-10)) +
           coord_flip() +
           annotation_custom(grobTree(textGrob(paste("speed =",format(as.numeric(1/sink_speed_regression$coefficients[1]), scientific = F, digits = 2)), x=0.68,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) +
           annotation_custom(grobTree(textGrob(expression(paste("  [",m~d^{-1},"]")), x=0.84,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=6, fontface="italic")))) +
           theme_bw()

         }
         
         
 #############################################################################
 ############################ RSUB PLOT ######################################
 #############################################################################
 ############################################################################# 
         
         DOXY_profile_data <- plotData_DOXY %>%
           filter(PARAM == "DOXY", cluster == i) %>%
           left_join(zp_data %>% select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
           left_join(mld_data %>% select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
           left_join(plotData_SIG %>% select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
           mutate(depths_below_zp = NA)
         DOXY_profile_data$zp <- round(DOXY_profile_data$zp)
         
         data_isodepth <- filter(DOXY_profile_data, jDay == 1) %>%
           filter(PRES >= zp)
         data_isodepth$depths_below_zp <- c(0:(1000-unique(data_isodepth$zp)))
         
         for (j in 2:365){
           data_by_day <- filter(DOXY_profile_data, jDay == j)
           data_isodepth_temp <- filter(data_by_day, PRES>=zp) 
           data_isodepth_temp$depths_below_zp <- c(0:(1000-unique(data_by_day$zp)))
           data_isodepth <- bind_rows(data_isodepth, data_isodepth_temp)
         }
         
  ########debut au climax
         # max_mld <- unique(DOXY_profile_data$jDay[which(DOXY_profile_data$mld==max(DOXY_profile_data$mld))]) +1
         # if(i==1){
         #   max_mld <- 82
         # }
         # if(i==2){
         #   max_mld <- 54
         # }
         # if(i==3){
         #   max_mld <- 56
         # }
         # if(i==4){
         #   max_mld <- 72
         # }
         # if(i==5){
         #   max_mld <- 76
         # }
         
         
    ########début à chaque max de O2  
         for (k in 1:(length(unique(data_isodepth$depths_below_zp)))){
            # print(k-1)
           data_isodepth_0_smth <- filter(data_isodepth, depths_below_zp==k-1)

           if(length(data_isodepth_0_smth$value) > 30){
            data_isodepth_0_smth <- filter(data_isodepth, depths_below_zp==k-1) %>%
               mutate(smth_values = smth(value, window = window, alpha = 2.5, method = "gaussian", tails = T))
            max_mld[k] <- match(max(data_isodepth_0_smth$smth_values, na.rm = T), data_isodepth_0_smth$smth_values)

            my_slope[k,] <- data_isodepth_0_smth %>%
              filter(jDay >= max_mld[k]) %>%
              select(dates = jDay, O2 = value, PRES = PRES, depths_below_zp = depths_below_zp) %>%
              group_by(depths_below_zp) %>% do(model = lm(O2~ dates, data = .)) %>% ungroup() %>% #groupe les données par "lvl" et applique un lm pour chaque "lvl". La fonction sors une liste pour chaque "lvl"
              transmute(depths_below_zp, coef = map(model, tidy)) %>% #la fonction tidy de broom permet de récupérer les infos de la list comme le coef std.error etc... Plour le R² modifier par la fonction "glance"
              unnest(coef) %>%
              filter(term == "dates") %>% #la fonction précédente retourne tout en une seule colonne. Unnest la sépare en plusieur colonne
              #optionel : filtrer les lignes relatives à l'intercept pour ne travailler que sur Y
              select(depths_below_zp,term, estimate, std.error) #optionnel: sélectionner que les "lvl" et intercept associés
           }

           if(length(data_isodepth_0_smth$value) < 30){
             max_mld[k] <- match(max(data_isodepth_0_smth$value, na.rm = T), data_isodepth_0_smth$value)
             my_slope[k,] <- data_isodepth_0_smth %>%
               filter(jDay >= max_mld[k]) %>%
               select(dates = jDay, O2 = value, PRES = PRES, depths_below_zp = depths_below_zp) %>%
               group_by(depths_below_zp) %>% do(model = lm(O2~ dates, data = .)) %>% ungroup() %>% #groupe les données par "lvl" et applique un lm pour chaque "lvl". La fonction sors une liste pour chaque "lvl"
               transmute(depths_below_zp, coef = map(model, tidy)) %>% #la fonction tidy de broom permet de récupérer les infos de la list comme le coef std.error etc... Plour le R² modifier par la fonction "glance"
               unnest(coef) %>%
               filter(term == "dates") %>% #la fonction précédente retourne tout en une seule colonne. Unnest la sépare en plusieur colonne
               #optionel : filtrer les lignes relatives à l'intercept pour ne travailler que sur Y
               select(depths_below_zp,term, estimate, std.error) #optionnel: sélectionner que les "lvl" et intercept associés
           }

         }


         # data_isodepth_0 <- filter(data_isodepth, depths_below_zp==0)
         # max_mld <- match(max(data_isodepth_0$value, na.rm = T), data_isodepth_0$value)
         
         #  ######## début à l'apex
         # print(max_mld)
         
         # if(i==1){
         #   max_mld <- 189
         # }
         # if(i==2){
         #   max_mld <- 178
         # }
         # if(i==3){
         #   max_mld <- 117
         # }
         # if(i==4){
         #   max_mld <- 89
         # }
         # if(i==5){
         #   max_mld <- 149
         # }

         # my_slope <- data_isodepth %>%
         #   filter(jDay >= max_mld) %>%
         #   select(dates = jDay, O2 = value, PRES = PRES, depths_below_zp = depths_below_zp) %>%
         #   group_by(depths_below_zp) %>% do(model = lm(O2~ dates, data = .)) %>% ungroup() %>% #groupe les données par "lvl" et applique un lm pour chaque "lvl". La fonction sors une liste pour chaque "lvl"
         #   transmute(depths_below_zp, coef = map(model, tidy)) %>% #la fonction tidy de broom permet de récupérer les infos de la list comme le coef std.error etc... Plour le R² modifier par la fonction "glance"
         #   unnest(coef) %>%
         #   filter(term == "dates") %>% #la fonction précédente retourne tout en une seule colonne. Unnest la sépare en plusieur colonne
         #   #optionel : filtrer les lignes relatives à l'intercept pour ne travailler que sur Y
         #   select(depths_below_zp,term, estimate, std.error) #optionnel: sélectionner que les "lvl" et intercept associés
         
         sigma <- mean(DOXY_profile_data$SIG) + 1000
         
         my_slope$estimate <- -(1/1.45)*my_slope$estimate *10^-3 *12 *sigma
         
         my_slope <-  my_slope %>% 
           mutate(sd = (1/1.45)*std.error *10^-3 *12 *sigma)
         
         R_plots[[i]] <- ggplot(data = my_slope, aes(x = depths_below_zp, y = estimate)) +
           geom_errorbar(aes(ymin = estimate-sd, ymax = estimate+sd), width = 15, color = "grey") +
           geom_line(size = 0.5, color = "black") +
           scale_x_reverse(name = "Pressure under Zp [dbar]", expand = c(0,0), limits = c(1000,-10)) +
           scale_y_continuous(name = expression(paste(net~respiration~"  [",mg~C~m^{-3}~d^{-1},"] ")), expand = c(0.02,0.02), limits = c(0, 1800)) +
           geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
           coord_flip() +
           theme_bw()
         
         
         resp_0_1000 <- select(my_slope, estimate) 
         resp_0_1000 <- as.vector(resp_0_1000$estimate)
         depths_below_zp_0_1000 <- select(my_slope, depths_below_zp)
         depths_below_zp_0_1000 <- as.vector(depths_below_zp_0_1000$depths_below_zp)
         sd_0_1000 <- select(my_slope, sd)
         sd_0_1000 <- as.vector(sd_0_1000$sd)
         
         resp_150_1000 <- filter(my_slope, depths_below_zp >= 150) %>%
           select(estimate)
         resp_150_1000 <- as.vector(resp_150_1000$estimate)
         depths_below_zp_150_1000 <- filter(my_slope, depths_below_zp >= 150) %>%
           select(depths_below_zp)
         depths_below_zp_150_1000 <- as.vector(depths_below_zp_150_1000$depths_below_zp)
         sd_150_1000 <- select(my_slope, sd)
         sd_150_1000 <- as.vector(sd_150_1000$sd)
         
         resp_300_1000 <- filter(my_slope, depths_below_zp >= 300) %>%
           select(estimate)
         resp_300_1000 <- as.vector(resp_300_1000$estimate)
         depths_below_zp_300_1000 <- filter(my_slope, depths_below_zp >= 300) %>%
           select(depths_below_zp)
         depths_below_zp_300_1000 <- as.vector(depths_below_zp_300_1000$depths_below_zp)
         sd_300_1000 <- select(my_slope, sd)
         sd_300_1000 <- as.vector(sd_300_1000$sd)
         
         resp_450_1000 <- filter(my_slope, depths_below_zp >= 450) %>%
           select(estimate)
         resp_450_1000 <- as.vector(resp_450_1000$estimate)
         depths_below_zp_450_1000 <- filter(my_slope, depths_below_zp >= 450) %>%
           select(depths_below_zp)
         depths_below_zp_450_1000 <- as.vector(depths_below_zp_450_1000$depths_below_zp)
         sd_450_1000 <- select(my_slope, sd)
         sd_450_1000 <- as.vector(sd_450_1000$sd)
         
         resp_600_1000 <- filter(my_slope, depths_below_zp >= 600) %>%
           select(estimate)
         resp_600_1000 <- as.vector(resp_600_1000$estimate)
         depths_below_zp_600_1000 <- filter(my_slope, depths_below_zp >= 600) %>%
           select(depths_below_zp)
         depths_below_zp_600_1000 <- as.vector(depths_below_zp_600_1000$depths_below_zp)
         sd_600_1000 <- select(my_slope, sd)
         sd_600_1000 <- as.vector(sd_600_1000$sd)
         
         resp_750_1000 <- filter(my_slope, depths_below_zp >= 750) %>%
           select(estimate)
         resp_750_1000 <- as.vector(resp_750_1000$estimate)
         depths_below_zp_750_1000 <- filter(my_slope, depths_below_zp >= 750) %>%
           select(depths_below_zp)
         depths_below_zp_750_1000 <- as.vector(depths_below_zp_750_1000$depths_below_zp)
         sd_750_1000 <- select(my_slope, sd)
         sd_750_1000 <- as.vector(sd_750_1000$sd)
         
         Rsub_0_1000 <- abs(depths_below_zp_0_1000[2]-depths_below_zp_0_1000[1])*abs(resp_0_1000[2]+resp_0_1000[1])/2
         sd_tot_0_1000 <- abs(depths_below_zp_0_1000[2]-depths_below_zp_0_1000[1])*abs(sd_0_1000[2]+sd_0_1000[1])/2
         for (k in 2:(length(my_slope$estimate)-1)){
           Rsub_0_1000 <- Rsub_0_1000 + abs(depths_below_zp_0_1000[k+1]-depths_below_zp_0_1000[k])*abs(resp_0_1000[k+1]+resp_0_1000[k])/2
           sd_tot_0_1000 <- sd_tot_0_1000 + abs(depths_below_zp_0_1000[k+1]-depths_below_zp_0_1000[k])*abs(sd_0_1000[k+1]+sd_0_1000[k])/2
         }
         
         Rsub_150_1000 <- abs(depths_below_zp_150_1000[2]-depths_below_zp_150_1000[1])*abs(resp_150_1000[2]+resp_150_1000[1])/2
         sd_tot_150_1000 <- abs(depths_below_zp_150_1000[2]-depths_below_zp_150_1000[1])*abs(sd_150_1000[2]+sd_150_1000[1])/2
         for (k in 2:(length(depths_below_zp_150_1000)-1)){
           Rsub_150_1000 <- Rsub_150_1000 + abs(depths_below_zp_150_1000[k+1]-depths_below_zp_150_1000[k])*abs(resp_150_1000[k+1]+resp_150_1000[k])/2
           sd_tot_150_1000 <- sd_tot_150_1000 + abs(depths_below_zp_150_1000[k+1]-depths_below_zp_150_1000[k])*abs(sd_150_1000[k+1]+sd_150_1000[k])/2
         }
         
         Rsub_300_1000 <- abs(depths_below_zp_300_1000[2]-depths_below_zp_300_1000[1])*abs(resp_300_1000[2]+resp_300_1000[1])/2
         sd_tot_300_1000 <- abs(depths_below_zp_300_1000[2]-depths_below_zp_300_1000[1])*abs(sd_300_1000[2]+sd_300_1000[1])/2
         for (k in 2:(length(depths_below_zp_300_1000)-1)){
           Rsub_300_1000 <- Rsub_300_1000 + abs(depths_below_zp_300_1000[k+1]-depths_below_zp_300_1000[k])*abs(resp_300_1000[k+1]+resp_300_1000[k])/2
           sd_tot_300_1000 <- sd_tot_300_1000 + abs(depths_below_zp_300_1000[k+1]-depths_below_zp_300_1000[k])*abs(sd_300_1000[k+1]+sd_300_1000[k])/2
         }
         
         Rsub_450_1000 <- abs(depths_below_zp_450_1000[2]-depths_below_zp_450_1000[1])*abs(resp_450_1000[2]+resp_450_1000[1])/2
         sd_tot_450_1000 <- abs(depths_below_zp_450_1000[2]-depths_below_zp_450_1000[1])*abs(sd_450_1000[2]+sd_450_1000[1])/2
         for (k in 2:(length(depths_below_zp_450_1000)-1)){
           Rsub_450_1000 <- Rsub_450_1000 + abs(depths_below_zp_450_1000[k+1]-depths_below_zp_450_1000[k])*abs(resp_450_1000[k+1]+resp_450_1000[k])/2
           sd_tot_450_1000 <- sd_tot_450_1000 + abs(depths_below_zp_450_1000[k+1]-depths_below_zp_450_1000[k])*abs(sd_450_1000[k+1]+sd_450_1000[k])/2
         }
         
         Rsub_600_1000 <- abs(depths_below_zp_600_1000[2]-depths_below_zp_600_1000[1])*abs(resp_600_1000[2]+resp_600_1000[1])/2
         sd_tot_600_1000 <- abs(depths_below_zp_600_1000[2]-depths_below_zp_600_1000[1])*abs(sd_600_1000[2]+sd_600_1000[1])/2
         for (k in 2:(length(depths_below_zp_600_1000)-1)){
           Rsub_600_1000 <- Rsub_600_1000 + abs(depths_below_zp_600_1000[k+1]-depths_below_zp_600_1000[k])*abs(resp_600_1000[k+1]+resp_600_1000[k])/2
           sd_tot_600_1000 <- sd_tot_600_1000 + abs(depths_below_zp_600_1000[k+1]-depths_below_zp_600_1000[k])*abs(sd_600_1000[k+1]+sd_600_1000[k])/2
         }
         
         Rsub_750_1000 <- abs(depths_below_zp_750_1000[2]-depths_below_zp_750_1000[1])*abs(resp_750_1000[2]+resp_750_1000[1])/2
         sd_tot_750_1000 <- abs(depths_below_zp_750_1000[2]-depths_below_zp_750_1000[1])*abs(sd_750_1000[2]+sd_750_1000[1])/2
         for (k in 2:(length(depths_below_zp_750_1000)-1)){
           Rsub_750_1000 <- Rsub_750_1000 + abs(depths_below_zp_750_1000[k+1]-depths_below_zp_750_1000[k])*abs(resp_750_1000[k+1]+resp_750_1000[k])/2
           sd_tot_750_1000 <- sd_tot_750_1000 + abs(depths_below_zp_750_1000[k+1]-depths_below_zp_750_1000[k])*abs(sd_750_1000[k+1]+sd_750_1000[k])/2
         }
         
         Rsub <- data.frame(Rsub = c(Rsub_0_1000, Rsub_150_1000, Rsub_300_1000, Rsub_450_1000, Rsub_600_1000, Rsub_750_1000)) %>%
           mutate(sd_tot = c(sd_tot_0_1000, sd_tot_150_1000, sd_tot_300_1000, sd_tot_450_1000, sd_tot_600_1000, sd_tot_750_1000)) %>%
           mutate(depths = c(0,150,300,450,600,750)) %>%
           mutate(cluster = rep(i,6))
         
         Rsub_plots[[i]] <- ggplot(data = Rsub, aes(x = depths, y = Rsub)) +
           #geom_point() +
           geom_line() +
           geom_errorbar(aes(ymin = Rsub-sd_tot, ymax = Rsub+sd_tot), width = 15) +
           scale_x_reverse(name = "layers below Zp [m] (y:1000m)", expand = c(0,0), limits = c(800,-10)) +
           scale_y_continuous(name = expression(paste(Integrated~net~respiration~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0.01,0.01), limits = c(0,1800)) +
           geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
           coord_flip() +
           theme_bw()
         
         begin_date_O2 <- c(mean(max_mld[1:150]),mean(max_mld[151:300]),mean(max_mld[301:450]),mean(max_mld[451:600]),mean(max_mld[601:750]))
         
         I_plots[[i]] <- I_plots[[i]] +
           
           geom_point(x = unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[1]))]), 
                      y = I_values_per_layers_ready_to_plot$depth_0_1000[which(I_values_per_layers_ready_to_plot$month == unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[1]))]))], color = "red", size = 0.5) +

           geom_point(x = unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[2]))]), 
                      y = I_values_per_layers_ready_to_plot$depth_150_1000[which(I_values_per_layers_ready_to_plot$month == unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[2]))]))], color = "red", size = 0.5)+
           
           geom_point(x = unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[3]))]), 
                      y = I_values_per_layers_ready_to_plot$depth_300_1000[which(I_values_per_layers_ready_to_plot$month == unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[3]))]))], color = "red", size = 0.5)+
           
           geom_point(x = unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[4]))]), 
                      y = I_values_per_layers_ready_to_plot$depth_450_1000[which(I_values_per_layers_ready_to_plot$month == unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[4]))]))], color = "red", size = 0.5)+
           
           geom_point(x = unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[5]))]), 
                      y = I_values_per_layers_ready_to_plot$depth_600_1000[which(I_values_per_layers_ready_to_plot$month == unique(I_values_per_layers_ready_to_plot$month[which(I_values_per_layers_ready_to_plot$dates == round(begin_date_O2[5]))]))], color = "red", size = 0.5)
         
 #############################################################################
 ############################ COMPARAISON ######################################
 #############################################################################
 #############################################################################        
        
        if(i==1){
          Rsub_Att <- bind_cols(att=slopes$slope, Rsub)
          regression_Rsub_Att[[i]] <- lm(data = Rsub_Att, formula = Rsub~att)
        }
        
        if(i!=1){
          Rsub_Att_temp <- bind_cols(att=slopes$slope, Rsub)
          regression_Rsub_Att[[i]] <- lm(data = Rsub_Att_temp, formula = Rsub~att)
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
        print(paste("cluster",i,": y = ",slope = regression_Rsub_Att[[i]]$coefficients[[2]],"x +",regression_Rsub_Att[[i]]$coefficients[[1]]))

   }
   
   Rsub_Att$cluster <- as.numeric(Rsub_Att$cluster)
   
   Rsub_Att_plots <- ggplot(data = Rsub_Att, aes(x = att, y = Rsub, group=cluster, fill=factor(cluster))) +
     geom_point(shape=21, size = 3) +
     scale_x_continuous(name = expression(paste(E[0]~" [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0.05,0.05)) +
     #scale_x_continuous(name = expression(paste("Transfert efficiency [%]")), expand = c(0.05,0.05)) +
     scale_y_continuous(name = expression(paste(R[sub]~~"[",mg~C~m^{-2}~d^{-1},"]")), expand = c(0.01,0.01)) +

     geom_abline(intercept = regression_Rsub_Att[[1]]$coefficients[[1]], slope = regression_Rsub_Att[[1]]$coefficients[[2]], color = "black", linetype = "dashed") +
     annotation_custom(grobTree(textGrob(expression(R^2~"="), x=0.70,  y=0.2, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[1]])$r.squared, scientific=F, digits = 3), x=0.79,  y=0.2, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     
     geom_abline(intercept = regression_Rsub_Att[[2]]$coefficients[[1]], slope = regression_Rsub_Att[[2]]$coefficients[[2]], color = "black", linetype = "dashed") +
     annotation_custom(grobTree(textGrob(expression(R^2~"="), x=0.70,  y=0.16, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[2]])$r.squared, scientific=F, digits = 3), x=0.79,  y=0.16, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     
     geom_abline(intercept = regression_Rsub_Att[[3]]$coefficients[[1]], slope = regression_Rsub_Att[[3]]$coefficients[[2]], color = "black", linetype = "dashed") +
     annotation_custom(grobTree(textGrob(expression(R^2~"="), x=0.70,  y=0.12, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[3]])$r.squared, scientific=F, digits = 3), x=0.79,  y=0.12, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     
     geom_abline(intercept = regression_Rsub_Att[[4]]$coefficients[[1]], slope = regression_Rsub_Att[[4]]$coefficients[[2]], color = "black", linetype = "dashed") +
     annotation_custom(grobTree(textGrob(expression(R^2~"="), x=0.70,  y=0.08, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[4]])$r.squared, scientific=F, digits = 3), x=0.79,  y=0.08, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     
     geom_abline(intercept = regression_Rsub_Att[[5]]$coefficients[[1]], slope = regression_Rsub_Att[[5]]$coefficients[[2]], color = "black", linetype = "dashed") +
     annotation_custom(grobTree(textGrob(expression(R^2~"="), x=0.70,  y=0.04, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     annotation_custom(grobTree(textGrob(format(summary(regression_Rsub_Att[[5]])$r.squared, scientific=F, digits = 3), x=0.79,  y=0.04, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
     
     theme_bw() +
     scale_fill_manual(name = "cluster",breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#FF6666","#FFCC00","#33CCFF", "#0066CC"),
                       labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) Austral_LCB","(4) PN","(5) AN")) 

    Rsub_Att <- Rsub_Att%>%
      mutate(ratio = Rsub / att) %>%
      mutate(clust_name = c(rep("(1) Austral_HCB",6),rep("(2) Austral_MCB",6),rep("(3) Austral_LCB",6),rep("(4) PN",6),rep("(5) AN",6))) %>%
      filter(att>0) 
    
   Rsub_Att$cluster <- as.factor(Rsub_Att$cluster)
   Rsub_Att$depths <- as.factor(Rsub_Att$depths)
   
   resp_plot <- ggplot(data = Rsub_Att, aes(x = clust_name, y = ratio, fill = clust_name)) +
     geom_boxplot(aes(fill = clust_name)) +
     theme_bw() +
     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, aes(fill = factor(depths))) +
     scale_fill_manual(values = c("#CC0000","#FF6666","#FFCC00","#33CCFF", "#0066CC",colors)) +
     theme(axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           legend.title = element_blank()) +
     labs(y="% of remineralization")
   
 layout <- c(
   area(1,1),
   area(1,2),
   area(1,3),
   area(2,1),
   area(2,2)
 )
 
 I_plots[[1]] <- I_plots[[1]] + theme(axis.title.x = element_blank())
 I_plots[[2]] <- I_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank())
 I_plots[[3]] <- I_plots[[3]] + theme(axis.title.y.left = element_blank())
 #I_plots[[4]] <- I_plots[[4]] 
 I_plots[[5]] <- I_plots[[5]] + theme(axis.title.y.left = element_blank())
 
 Ez_plots[[1]] <- Ez_plots[[1]] + theme(axis.title.x = element_blank())
 Ez_plots[[2]] <- Ez_plots[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
 Ez_plots[[3]] <- Ez_plots[[3]] + theme(axis.title.y = element_blank())
 #Ez_plots[[4]] <- Ez_plots[[4]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
 Ez_plots[[5]] <- Ez_plots[[5]] + theme(axis.title.y = element_blank())
 
 Rp_plots[[1]] <- Rp_plots[[1]] + theme(axis.title.x = element_blank())
 Rp_plots[[2]] <- Rp_plots[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
 Rp_plots[[3]] <- Rp_plots[[3]] + theme(axis.title.y = element_blank())
 #Rp_plots[[4]] <- Rp_plots[[4]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
 Rp_plots[[5]] <- Rp_plots[[5]] + theme(axis.title.y = element_blank())
 
 Att_plots[[1]] <- Att_plots[[1]] + theme(axis.title.x = element_blank())
 Att_plots[[2]] <- Att_plots[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
 Att_plots[[3]] <- Att_plots[[3]] + theme(axis.title.y = element_blank())
 #Att_plots[[4]] <- Att_plots[[4]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
 Att_plots[[5]] <- Att_plots[[5]] + theme(axis.title.y = element_blank())
 
 Instant_flux_plots[[1]] <- Instant_flux_plots[[1]] + theme(axis.title.x = element_blank())
 Instant_flux_plots[[2]] <- Instant_flux_plots[[2]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
 Instant_flux_plots[[3]] <- Instant_flux_plots[[3]] + theme(axis.title.y = element_blank())
 #Instant_flux_plots[[4]] <- Instant_flux_plots[[4]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
 Instant_flux_plots[[5]] <- Instant_flux_plots[[5]] + theme(axis.title.y = element_blank())
 
 # Instant_norm_flux_plots[[1]] <- Instant_norm_flux_plots[[1]] + theme(axis.title.x = element_blank())
 # Instant_norm_flux_plots[[2]] <- Instant_norm_flux_plots[[2]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank())
 # Instant_norm_flux_plots[[3]] <- Instant_norm_flux_plots[[3]] + theme(axis.title.y = element_blank())
 # #Instant_norm_flux_plots[[4]] <- Instant_norm_flux_plots[[4]] + theme(axis.title.y = element_blank(),axis.title.x = element_blank())
 # Instant_norm_flux_plots[[5]] <- Instant_norm_flux_plots[[5]] + theme(axis.title.y = element_blank())
 
 Rsub_plots[[1]] <- Rsub_plots[[1]] + theme(axis.title.x = element_blank())
 Rsub_plots[[2]] <- Rsub_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank())
 Rsub_plots[[3]] <- Rsub_plots[[3]] + theme(axis.title.y.left = element_blank())
 #Rsub_plots[[4]] <- Rsub_plots[[4]] 
 Rsub_plots[[5]] <- Rsub_plots[[5]] + theme(axis.title.y.left = element_blank())
 
 # Rsub_Att_plots[[1]] <- Rsub_Att_plots[[1]] + theme(axis.title.x = element_blank())
 # Rsub_Att_plots[[2]] <- Rsub_Att_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank())
 # Rsub_Att_plots[[3]] <- Rsub_Att_plots[[3]] + theme(axis.title.y.left = element_blank())
 # #Rsub_Att_plots[[4]] <- Rsub_Att_plots[[4]]
 # Rsub_plots[[5]] <- Rsub_Att_plots[[5]] + theme(axis.title.y.left = element_blank())
 
 finalPlot_I <- wrap_plots(I_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 finalPlot_Ez <- wrap_plots(Ez_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 finalPlot_Rp <- wrap_plots(Rp_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 finalPlot_Att <- wrap_plots(Att_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 finalPlot_Instant <- wrap_plots(Instant_flux_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 # finalPlot_norm_Instant <- wrap_plots(Instant_norm_flux_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 finalPlot_Rsub <- wrap_plots(Rsub_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 # finalPlot_Rsub_Att <- wrap_plots(Rsub_Att_plots,design=layout) #+ plot_annotation(tag_levels = "1")
 
 finalPlot_Sink <- sinking_speed_plots[[1]] + sinking_speed_plots[[5]]
 #finalPlot_Sink <- wrap_plots(sinking_speed_plots,ncol = 2) #+ plot_annotation(tag_levels = "1")
 
 ggsave(plot = finalPlot_I,"figures/I_plots.png", width = 11.5, height = 5, dpi = 300)
 ggsave(plot = finalPlot_Ez,"figures/Ez_plots.png", width = 11.5, height = 5, dpi = 300)
 ggsave(plot = finalPlot_Rp,"figures/Rp_plots.png", width = 11.5, height = 5, dpi = 300)
 ggsave(plot = finalPlot_Att,"figures/Att_plots.png", width = 11, height = 8, dpi = 300)
 ggsave(plot = finalPlot_Instant,"figures/Instant_flux_plots.png", width = 11, height = 8, dpi = 300)
 # ggsave(plot = finalPlot_norm_Instant,"figures/Instant_flux_norm_plots.png", width = 11, height = 8, dpi = 300)
 ggsave(plot = finalPlot_Rsub,"figures/Rsub_plots.png", width = 11, height = 8, dpi = 300)
 ggsave(plot = finalPlot_Sink,"figures/Sink_plots.png", width = 6.2, height = 4, dpi = 300)
 # ggsave(plot = finalPlot_Rsub_Att,"figures/Rsub_Att_plots.png", width = 11, height = 8, dpi = 300)
 ggsave(plot = Rsub_Att_plots,"figures/Rsub_Att_plots.png", width = 6, height = 4, dpi = 300)
 ggsave(plot = resp_plot,"figures/ratio_resp.png", width = 6, height = 4, dpi = 300)
 } 
 
##############################################################################
 
 Kheireddine_2020_processing_AOU <- function(values, depths, dates, zp, depth_bin, mld, metrics, window = 30) {
   
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
   Att_plots <- list()
   names <- c("(1) Austral HCB","(2) Austral MCB","(3) STZ","(4) PN","(5) AN")
   
   for (i in 1:length(unique(profile_data_DOXY$cluster))){
      print(paste("cluster", i))
     
       AOU_profile_data <- profile_data_DOXY %>%
         filter(PARAM == "AOU", cluster == i) %>%
         left_join(zp_data %>% select(jDay, PRES, zp = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
         left_join(SIG_data %>% select(jDay, PRES, SIG = value, cluster), by = c("jDay", "PRES", "cluster")) %>%
         left_join(mld_data %>% select(jDay, PRES, mld = value, cluster), by = c("jDay", "PRES", "cluster"))
       
       #### Arguments
       
       depths = AOU_profile_data$PRES
       dates = AOU_profile_data$jDay
       zp = AOU_profile_data$zp
       mld = AOU_profile_data$mld
       depth_bin = 100
       SIG <- AOU_profile_data$SIG +1000
       values <- AOU_profile_data$value * SIG
       # metrics = read_csv(paste("data/metrics_cluster",i,".csv", sep=""), show_col_types = F)
     
       depth_seq <- c("whole_column", "productive_layer", seq(0, 800-depth_bin, by = depth_bin))
       # depth_seq_below_PPZ <- which( numbers_only(depth_seq) )
       # 
       depth_layers <- c("whole_column","productive_layer","depth_0_100","depth_100_200","depth_200_300","depth_300_400","depth_400_500","depth_500_600","depth_600_700","depth_700_800" )       # deepest_depth_layer <- depth_layers[which(depth_seq == max_depth_to_consider)]
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
           (mean(df$values[index]))
         })
         
         data.frame(df[1,], layer_name = depth_layers, DEPTH = depth_seq, I = integrated_values, row.names = NULL)
         
       }, .parallel = FALSE, .id = NULL) %>% 
         dplyr::group_by(layer_name) %>%
         dplyr::mutate(smoothed_I = smth(I, window = window, alpha = 2.5, method = "gaussian", tails = T),#*0.001,
                Ez_value = c( NA, diff(smoothed_I) / {diff(dates)})) %>% 
         ungroup()
       
       I_values_per_layers_ready_to_plot <- values_per_layers %>% pivot_wider(id_cols = c(dates, zp, mld), 
                                                                              names_from = "layer_name",values_from = "smoothed_I")
       
       I_plots[[i]] <- ggplot(data = I_values_per_layers_ready_to_plot/1000, aes(x = dates*1000)) 
       
       # geom_rect(data = flux_dates_df, aes(xmin = min_date, xmax = max_date, ymin = min_y, ymax = max_y), fill = "green3", alpha = 0.25)
       # geom_area(aes(y = whole_column), fill = "yellow", alpha = 0.5, color = "black", size = 0.15)
       
       depth_nb <- depth_seq[numbers_only(depth_seq)] %>% as.numeric()
       colors <- paste("gray", rescale(x =  depth_nb, to = c(0,80)) %>% round(), sep = "")
       
       for (j in rev(seq_along(depth_nb)) ) {
         I_plots[[i]] <- I_plots[[i]] + geom_area(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), fill = colors[j], color = colors[j]) +
           geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "black", size = 0.15) +
           geom_line(aes_string(y = grep(x = depth_layers, pattern = "depth", value = TRUE)[j]), color = "white", linetype = "dashed", size = 0.15)
       }
       
       
       # onset_x = metrics[1,1] %>% as.numeric()
       # climax_x = metrics[2,1] %>% as.numeric()
       # apex_x = metrics[3,1] %>% as.numeric()
       # min_0_1000_y = I_values_per_layers_ready_to_plot$depth_0_1000[apex_x]
       # max_0_1000_y = I_values_per_layers_ready_to_plot$depth_0_1000[365]
       # slope <- (max_0_1000_y-min_0_1000_y) / (365-apex_x)
       
       # trace bloom phenology with iPOC 0_1000
       whole_column <- values_per_layers %>%
         filter(DEPTH == "whole_column")
       iAOU_whole_column <- whole_column$smoothed_I
       
       depth_0_100 <- values_per_layers %>%
         filter(DEPTH == "0")
       iAOU_0_100 <- depth_0_100$smoothed_I
       
       depth_100_200 <- values_per_layers %>%
         filter(DEPTH == "100")
       iAOU_100_200 <- depth_100_200$smoothed_I
       
       depth_200_300 <- values_per_layers %>%
         filter(DEPTH == "200")
       iAOU_200_300 <- depth_200_300$smoothed_I
       
       depth_300_400 <- values_per_layers %>%
         filter(DEPTH == "300")
       iAOU_300_400 <- depth_300_400$smoothed_I
       
       depth_400_500 <- values_per_layers %>%
         filter(DEPTH == "400")
       iAOU_400_500 <- depth_400_500$smoothed_I
       
       depth_500_600 <- values_per_layers %>%
         filter(DEPTH == "500")
       iAOU_500_600 <- depth_500_600$smoothed_I
       
       depth_600_700 <- values_per_layers %>%
         filter(DEPTH == "600")
       iAOU_600_700 <- depth_600_700$smoothed_I
       
       depth_700_800 <- values_per_layers %>%
         filter(DEPTH == "700")
       iAOU_700_800 <- depth_700_800$smoothed_I
       
      
           ######### 0
           min_0_100_x <- match(min(iAOU_0_100),iAOU_0_100)
           min_0_100_y <- iAOU_0_100[min_0_100_x]
           max_0_100_x <- match(max(iAOU_0_100),iAOU_0_100)
           max_0_100_y <- iAOU_0_100[max_0_100_x]
           slope_0_100 <- (max_0_100_y - min_0_100_y) / (max_0_100_x + 365 - min_0_100_x) 
           ######### 150
           min_100_200_x <- match(min(iAOU_100_200),iAOU_100_200)
           min_100_200_y <- iAOU_100_200[min_100_200_x]
           max_100_200_x <- match(max(iAOU_100_200),iAOU_100_200)
           max_100_200_y <- iAOU_100_200[max_100_200_x]
           slope_100_200 <- (max_100_200_y - min_100_200_y) / (max_100_200_x + 365 - min_100_200_x)
           ######### 300
           min_200_300_x <- match(min(iAOU_200_300),iAOU_200_300)
           min_200_300_y <- iAOU_200_300[min_200_300_x]
           max_200_300_x <- match(max(iAOU_200_300),iAOU_200_300)
           max_200_300_y <- iAOU_200_300[max_200_300_x]
           slope_200_300 <- (max_200_300_y - min_200_300_y) / (max_200_300_x + 365 - min_200_300_x)
           ######### 450
           min_300_400_x <- match(min(iAOU_300_400),iAOU_300_400)
           min_300_400_y <- iAOU_300_400[min_300_400_x]
           max_300_400_x <- match(max(iAOU_300_400),iAOU_300_400)
           max_300_400_y <- iAOU_300_400[max_300_400_x]
           slope_300_400 <- (max_300_400_y - min_300_400_y) / (max_300_400_x + 365 - min_300_400_x)
           ######### 600
           min_400_500_x <- match(min(iAOU_400_500),iAOU_400_500)
           min_400_500_y <- iAOU_400_500[min_400_500_x]
           max_400_500_x <- match(max(iAOU_400_500),iAOU_400_500)
           max_400_500_y <- iAOU_400_500[max_400_500_x]
           slope_400_500 <- (max_400_500_y - min_400_500_y) / (max_400_500_x + 365 - min_400_500_x)
           ######### 750
           min_500_600_x <- match(min(iAOU_500_600),iAOU_500_600)
           min_500_600_y <- iAOU_500_600[min_500_600_x]
           max_500_600_x <- match(max(iAOU_500_600),iAOU_500_600)
           max_500_600_y <- iAOU_500_600[max_500_600_x]
           slope_500_600 <- (max_500_600_y - min_500_600_y) / (max_500_600_x + 365 - min_500_600_x)
           
           min_600_700_x <- match(min(iAOU_600_700),iAOU_600_700)
           min_600_700_y <- iAOU_600_700[min_600_700_x]
           max_600_700_x <- match(max(iAOU_600_700),iAOU_600_700)
           max_600_700_y <- iAOU_600_700[max_600_700_x]
           slope_600_700 <- (max_600_700_y - min_600_700_y) / (max_600_700_x + 365 - min_600_700_x)
           
           min_700_800_x <- match(min(iAOU_700_800),iAOU_700_800)
           min_700_800_y <- iAOU_700_800[min_700_800_x]
           max_700_800_x <- match(max(iAOU_700_800),iAOU_700_800)
           max_700_800_y <- iAOU_700_800[max_700_800_x]
           slope_700_800 <- (max_700_800_y - min_700_800_y) / (max_700_800_x + 365 - min_700_800_x)
       
       
       slopes <- data.frame(flux = c(slope_0_100,slope_100_200,
                                     slope_200_300,slope_300_400,
                                     slope_400_500,slope_500_600,
                                     slope_600_700,slope_700_800))
       # %>%
       #   mutate(i_depths = depth_nb+1)
       slopes$flux[1] <- -(1/1.45)*slopes$flux[1] *10^-3 *12 *100  # x le quotient respiratoire pour passer de O2 en C, 
       slopes$flux[2] <- -(1/1.45)*slopes$flux[2] *10^-3 *12 *100  # x0.001 pour passer de µmol en mmol,
       slopes$flux[3] <- -(1/1.45)*slopes$flux[3] *10^-3 *12 *100  # x12 pour passer mmol de C en mg de C
       slopes$flux[4] <- -(1/1.45)*slopes$flux[4] *10^-3 *12 *100  # x100 pour integrer sur une couche de 100m
       slopes$flux[5] <- -(1/1.45)*slopes$flux[5] *10^-3 *12 *100 
       slopes$flux[6] <- -(1/1.45)*slopes$flux[6] *10^-3 *12 *100 
       
       if(i!=5){
         slopes$flux[7] <- -(1/1.45)*slopes$flux[7] *10^-3 *12 *100
         slopes$flux[8] <- -(1/1.45)*slopes$flux[8] *10^-3 *12 *100
       }
       
       slopes <- slopes %>%
         mutate(prof = c(50,150,250,350,450,550,650,750))
       
       if(i==5){
         slopes <- slopes %>% filter(prof < 650)
         # incert <- filter(incertitudes_C, cluster==i) %>%
         #   mutate(prof = c(50,150,250,350,450,550))
       }
       # 
       # if(i!=5){
       #   
       #   # incert <- filter(incertitudes_C, cluster==i) %>%
       #   #   mutate(prof = c(50,150,250,350,450,550,650,750))
       # }
       
       scaling_factor_3 <- sec_axis_adjustement_factors(-c(0,600), c(-20,300))
       
       
       I_plots[[i]] <- print(I_plots[[i]] +
         # geom_line(aes(y = productive_layer, color = "y1"), size=0.5) + 
         scale_x_continuous(name = "relative month", expand = c(0,0)) +
         scale_y_continuous(expression(paste(Average~AOU,"  [",µmol~kg^{-1},"]")), breaks =  c(-20, 0, 50, 100, 150, 200, 250, 300), sec.axis = sec_axis(trans = ~ {. - scaling_factor_3$adjust} / scaling_factor_3$diff *-1, name = "MLD [m]")) + 
         geom_line(aes(y = -I_values_per_layers_ready_to_plot$mld * scaling_factor_3$diff + scaling_factor_3$adjust, color = "y2"), size = 0.5) + 
         scale_color_manual(values = c("y1" = "#CC0000", "y2" = "#3399FF")) +
         coord_cartesian(xlim = c(0,365), ylim = c(-20,300)) + 
         #geom_point(aes(x = apex_x, y = apex_y), shape=21, fill = "#FF3300", size=2, color="#FF3300") +
         # geom_vline(xintercept = apex_x_month, linetype = "dotted", color = "#FF0000", size = 1.3) +
         # geom_point(aes(x = onset_x, y = onset_y), shape=21, fill = "#FFCC33", size=2, color="#FFCC33") +
         # geom_vline(xintercept = onset_x_month, linetype = "dotted", color = "#339900", size = 1.3) +
         # geom_point(aes(x = climax_x, y = climax_y), shape=21, fill = "#FF9900", size=2, color="#FF9900") +
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
         ggtitle(names[i]))
       
       ##########################################################################################################################################################################################################
       ##########################################################################################################################################################################################################
       ##########################################################################################################################################################################################################
       ##########################################################################################################################################################################################################
       ##########################################################################################################################################################################################################
       
       if(i!=5){
         Att_plots[[i]] <- print(ggplot(data = slopes, aes(x = prof, y = -flux)) +
                                   geom_col(col="black",fill = colors) +
                                   # geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
                                   # geom_point(shape=21, fill = colors, size=2) +
                                   # geom_function(fun = f, args = list(a=a, b=b), color = "darkred", size = 1) +
                                   # geom_line(data = regression_data, aes(x = z, y = flux), color ="darkred", size = 1) +
                                   scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(700,0)) +
                                   scale_y_continuous(name = expression(paste(Seasonal~accumulation,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,700)) +
                                   coord_flip() +
                                   # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                                   theme_bw() +
                                   theme(axis.text = element_text(size = 11, colour = "black"),
                                         axis.title = element_text(size = 11, colour = "black"),
                                         axis.text.x=element_text(angle=0),
                                         plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                                   # annotation_custom(grobTree(textGrob(paste("b = ",format(b, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   
                                   ggtitle(names[i]))
       }
       
       if(i==5){
         Att_plots[[i]] <- print(ggplot(data = slopes, aes(x = prof, y = -flux)) +
                                   geom_col(col="black",fill = colors[-(7:8)]) +
                                   # geom_errorbar(aes(ymin = slope_min, ymax = slope_max), width = 18, col="darkgrey") +
                                   # geom_point(shape=21, fill = colors[-(7:8)], size=2) +
                                   # xlim(c(0,40)) +
                                   # geom_function(fun = f, args = list(a=a, b=b), color = "darkred",size = 1) +
                                   # geom_line(data = regression_data[(1:551),], aes(x = z, y = flux), color ="darkred", size = 1) +
                                   scale_x_reverse(name = "profondeur en dessous de Zp [m]", expand = c(0,0), limits = c(700,0)) +
                                   scale_y_continuous(name = expression(paste(Accumulation~saisonnière,"  [",mg~C~m^{-2}~d^{-1},"]")), expand = c(0,0),limits = c(0,700)) +
                                   coord_flip() +
                                   # geom_line(data = regression_data, aes(x=z, y=pred_slope), color = 'red')+
                                   theme_bw() +
                                   theme(axis.text = element_text(size = 11, colour = "black"),
                                         axis.title = element_text(size = 11, colour = "black"),
                                         axis.text.x=element_text(angle=0),
                                         plot.title = element_text(size=13, face="bold.italic", hjust = 0.5)) +
                                   # annotation_custom(grobTree(textGrob(paste("b = ",format(b, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("c = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.75, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("l = ",format(b, scientific = F, digits = 4)), x=0.56,  y=0.70, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("C = ",format(c, scientific = F, digits = 3)), x=0.56,  y=0.65, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("max = ",format(slopes$slope[1], scientific = F, digits = 3)), x=0.56,  y=0.60, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   # annotation_custom(grobTree(textGrob(paste("min = ",format(slopes$slope[6], scientific = F, digits = 3)), x=0.56,  y=0.55, hjust=0, gp=gpar(col="black", fontsize=10, fontface="italic")))) + 
                                   
                                   ggtitle(names[i]))
       }
       
   }
   
   return(list("I_plot" = I_plot,
               # "Ez_plot" = Ez_plot,
               data = values_per_layers))
   
 }
 
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
 

 #####################
   