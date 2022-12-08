R_plots <- list()
Rsub_plots <- list()
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


for (i in 1:5){

      print(paste("cluster ",i))
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
          # for (l in 1:round((1000-unique(data_by_day$zp)))){
          #   data_isodepth_temp <- filter(data_by_day, PRES == (zp + l))
          #   data_isodepth_temp$depths_below_zp <- l
          #   data_isodepth <- bind_rows(data_isodepth, data_isodepth_temp)
          #   }
      }
        
      max_mld <- unique(DOXY_profile_data$jDay[which(DOXY_profile_data$mld==max(DOXY_profile_data$mld))]) +1
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
        # 
      my_slope <- data_isodepth %>%
        filter(jDay >= max_mld) %>%
        select(dates = jDay, O2 = value, PRES = PRES, depths_below_zp = depths_below_zp) %>%
        group_by(depths_below_zp) %>% do(model = lm(O2~ dates, data = .)) %>% ungroup() %>% #groupe les données par "lvl" et applique un lm pour chaque "lvl". La fonction sors une liste pour chaque "lvl"
        transmute(depths_below_zp, coef = map(model, tidy)) %>% #la fonction tidy de broom permet de récupérer les infos de la list comme le coef std.error etc... Plour le R² modifier par la fonction "glance"
        unnest(coef) %>%
        filter(term == "dates") %>% #la fonction précédente retourne tout en une seule colonne. Unnest la sépare en plusieur colonne
        #optionel : filtrer les lignes relatives à l'intercept pour ne travailler que sur Y
        select(depths_below_zp,term, estimate, std.error) #optionnel: sélectionner que les "lvl" et intercept associés
      
      sigma <- mean(DOXY_profile_data$SIG) + 1000
      
      my_slope$estimate <- -(1/1.45)*my_slope$estimate *10^-3 *12 *sigma
      
      my_slope <-  my_slope %>% 
        mutate(sd = (1/1.45)*std.error *10^-3 *12 *sigma)
      
      R_plots[[i]] <- ggplot(data = my_slope, aes(x = depths_below_zp, y = estimate)) +
        geom_line(size = 1, color = "black") +
        geom_errorbar(aes(ymin = estimate-sd, ymax = estimate+sd), width = 15) +
        scale_x_reverse(name = "Pressure under Zp [dbar]", expand = c(0,0), limits = c(1000,-10)) +
        scale_y_continuous(name = expression(paste(net~respiration~"  [",mg~C~m^{-3}~d^{-1},"] ")), expand = c(0.02,0.02), limits = c(-0.15, 2.1)) +
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
      
      Rsub_0_1000[i] <- abs(depths_below_zp_0_1000[2]-depths_below_zp_0_1000[1])*abs(resp_0_1000[2]-resp_0_1000[1])
      sd_tot_0_1000[i] <- abs(depths_below_zp_0_1000[2]-depths_below_zp_0_1000[1])*abs(sd_0_1000[2]-sd_0_1000[1])
      for (k in 2:(length(my_slope$estimate)-1)){
        Rsub_0_1000[i] <- Rsub_0_1000[i] + abs(depths_below_zp_0_1000[k+1]-depths_below_zp_0_1000[k])*abs(resp_0_1000[k+1]-resp_0_1000[k])
        sd_tot_0_1000[i] <- sd_0_1000[i] + abs(depths_below_zp_0_1000[k+1]-depths_below_zp_0_1000[k])*abs(sd_0_1000[k+1]-sd_0_1000[k])
      }
      
      Rsub_150_1000[i] <- abs(depths_below_zp_150_1000[2]-depths_below_zp_150_1000[1])*abs(resp_150_1000[2]-resp_150_1000[1])
      sd_tot_150_1000[i] <- abs(depths_below_zp_150_1000[2]-depths_below_zp_150_1000[1])*abs(sd_150_1000[2]-sd_150_1000[1])
      for (k in 2:(length(depths_below_zp_150_1000)-1)){
        Rsub_150_1000[i] <- Rsub_150_1000[i] + abs(depths_below_zp_150_1000[k+1]-depths_below_zp_150_1000[k])*abs(resp_150_1000[k+1]-resp_150_1000[k])
        sd_tot_150_1000[i] <- sd_150_1000[i] + abs(depths_below_zp_150_1000[k+1]-depths_below_zp_150_1000[k])*abs(sd_150_1000[k+1]-sd_150_1000[k])
      }
      
      Rsub_300_1000[i] <- abs(depths_below_zp_300_1000[2]-depths_below_zp_300_1000[1])*abs(resp_300_1000[2]-resp_300_1000[1])
      sd_tot_300_1000[i] <- abs(depths_below_zp_300_1000[2]-depths_below_zp_300_1000[1])*abs(sd_300_1000[2]-sd_300_1000[1])
      for (k in 2:(length(depths_below_zp_300_1000)-1)){
        Rsub_300_1000[i] <- Rsub_300_1000[i] + abs(depths_below_zp_300_1000[k+1]-depths_below_zp_300_1000[k])*abs(resp_300_1000[k+1]-resp_300_1000[k])
        sd_tot_300_1000[i] <- sd_300_1000[i] + abs(depths_below_zp_300_1000[k+1]-depths_below_zp_300_1000[k])*abs(sd_300_1000[k+1]-sd_300_1000[k])
      }
      
      Rsub_450_1000[i] <- abs(depths_below_zp_450_1000[2]-depths_below_zp_450_1000[1])*abs(resp_450_1000[2]-resp_450_1000[1])
      sd_tot_450_1000[i] <- abs(depths_below_zp_450_1000[2]-depths_below_zp_450_1000[1])*abs(sd_450_1000[2]-sd_450_1000[1])
      for (k in 2:(length(depths_below_zp_450_1000)-1)){
        Rsub_450_1000[i] <- Rsub_450_1000[i] + abs(depths_below_zp_450_1000[k+1]-depths_below_zp_450_1000[k])*abs(resp_450_1000[k+1]-resp_450_1000[k])
        sd_tot_450_1000[i] <- sd_450_1000[i] + abs(depths_below_zp_450_1000[k+1]-depths_below_zp_450_1000[k])*abs(sd_450_1000[k+1]-sd_450_1000[k])
      }
      
      Rsub_600_1000[i] <- abs(depths_below_zp_600_1000[2]-depths_below_zp_600_1000[1])*abs(resp_600_1000[2]-resp_600_1000[1])
      sd_tot_600_1000[i] <- abs(depths_below_zp_600_1000[2]-depths_below_zp_600_1000[1])*abs(sd_600_1000[2]-sd_600_1000[1])
      for (k in 2:(length(depths_below_zp_600_1000)-1)){
        Rsub_600_1000[i] <- Rsub_600_1000[i] + abs(depths_below_zp_600_1000[k+1]-depths_below_zp_600_1000[k])*abs(resp_600_1000[k+1]-resp_600_1000[k])
        sd_tot_600_1000[i] <- sd_600_1000[i] + abs(depths_below_zp_600_1000[k+1]-depths_below_zp_600_1000[k])*abs(sd_600_1000[k+1]-sd_600_1000[k])
      }
      
      Rsub_750_1000[i] <- abs(depths_below_zp_750_1000[2]-depths_below_zp_750_1000[1])*abs(resp_750_1000[2]-resp_750_1000[1])
      sd_tot_750_1000[i] <- abs(depths_below_zp_750_1000[2]-depths_below_zp_750_1000[1])*abs(sd_750_1000[2]-sd_750_1000[1])
      for (k in 2:(length(depths_below_zp_750_1000)-1)){
        Rsub_750_1000[i] <- Rsub_750_1000[i] + abs(depths_below_zp_750_1000[k+1]-depths_below_zp_750_1000[k])*abs(resp_750_1000[k+1]-resp_750_1000[k])
        sd_tot_750_1000[i] <- sd_750_1000[i] + abs(depths_below_zp_750_1000[k+1]-depths_below_zp_750_1000[k])*abs(sd_750_1000[k+1]-sd_750_1000[k])
      }
      
}
      
Rsub <- data.frame(Rsub = c(Rsub_0_1000, Rsub_150_1000, Rsub_300_1000, Rsub_450_1000, Rsub_600_1000, Rsub_750_1000)) %>%
  mutate(sd_tot = c(sd_tot_0_1000, sd_tot_150_1000, sd_tot_300_1000, sd_tot_450_1000, sd_tot_600_1000, sd_tot_750_1000)) %>%
  mutate(depths = c(rep(0,5),rep(150,5), rep(300,5), rep(450,5), rep(600,5), rep(750,5))) %>%
  mutate(cluster = rep(c("Austral HCB", "Austral MCB", "Austral LCB", "PN", "AN"),6)) %>%
  mutate(cluster_nb = rep(c("1", "2", "3", "4", "5"),6)) 


for (i in 1:5){
  
  Rsub_temp <- filter(Rsub, cluster_nb == i)
  
Rsub_plots[[i]] <- ggplot(data = Rsub_temp, aes(x = depths, y = Rsub)) +
  #geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = Rsub-sd_tot, ymax = Rsub+sd_tot), width = 15) +
  scale_x_reverse(name = "Depths under zp [dbar]", expand = c(0,0), limits = c(800,-10)) +
  scale_y_continuous(name = expression(paste(Integrated~net~respiration~"  [",mg~C~m^{-2}~d^{-1},"] ")), expand = c(0.1,0.1), limits = c(0,3.3)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
  coord_flip() +
  theme_bw()

}

layout <- c(
  area(1,1),
  area(1,2),
  area(1,3),
  area(2,1),
  area(2,2)
)

Rsub_plots[[1]] <- Rsub_plots[[1]] + theme(axis.title.x = element_blank())
Rsub_plots[[2]] <- Rsub_plots[[2]] + theme(axis.title.y.left = element_blank(), axis.title.x = element_blank())
Rsub_plots[[3]] <- Rsub_plots[[3]] + theme(axis.title.y.left = element_blank())
#Rsub_plots[[4]] <- Rsub_plots[[4]] 
Rsub_plots[[5]] <- Rsub_plots[[5]] + theme(axis.title.y.left = element_blank())

finalPlot_Rsub <- wrap_plots(Rsub_plots,design=layout) #+ plot_annotation(tag_levels = "1")
ggsave(plot = finalPlot_Rsub,"figures/Rsub_plots.png", width = 11, height = 8, dpi = 300)


Austral_MCB <- filter(Rsub, cluster_nb==2)
plot(depth_1$jDay, depth_1$value)
