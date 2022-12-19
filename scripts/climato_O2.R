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


plotData_DOXY <- read_csv("data/DOXY/plotData.csv") 
svdMatrix <- read_csv("data/DOXY/svdMatrix_timeseries_DOXY.csv") %>% filter(PARAM=="DOXY")


climato <- plotData_DOXY %>%
  filter(PARAM == "DOXY") %>%
  group_by(cluster,PRES) %>%
  dplyr::summarise(cluster = unique(cluster),
                   PRES = unique(PRES), 
                   climato = mean(value),
                   sd = sd(value))


resp_each_depth_plot <- ggplot(data = climato, aes(x = PRES, y = climato, group = cluster, fill = factor(cluster), colour = factor(cluster))) +
  geom_line() +
  # geom_errorbar(aes(ymin = estimate-sd, ymax = estimate+sd), width = 15) +
  scale_x_reverse(name = "PRES", expand = c(0,0), limits = c(850,-10)) +
  scale_y_continuous(name = expression(paste(climato~O2~" [",µmol~O2~kg^{-1},"] ")), expand = c(0.02,0.02), limits = c(0, 350)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = 0.8) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                    labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN")) +
  scale_colour_manual(name = element_blank(), breaks = c("1", "2", "3","4","5"), values=c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#0066CC"),
                      labels = c("(1) Austral_HCB","(2) Austral_MCB","(3) STZ","(4) PN","(5) AN"))

inputMatrix <- inputMatrix %>% filter(PARAM== "DOXY")

climato_1 <- climato %>% 
  filter(cluster==1) %>% 
  select(climato)
climato_1 <- as.numeric(climato_1$climato)
climato_1 <- rep(climato_1,365)

climato_2 <- climato %>% 
  filter(cluster==2) %>% 
  select(climato)
climato_2 <- as.numeric(climato_2$climato)
climato_2 <- rep(climato_2,365)

climato_3 <- climato %>% 
  filter(cluster==3) %>% 
  select(climato)
climato_3 <- as.numeric(climato_3$climato)
climato_3 <- rep(climato_3,365)

climato_4 <- climato %>% 
  filter(cluster==4) %>% 
  select(climato)
climato_4 <- as.numeric(climato_4$climato)
climato_4 <- rep(climato_4,365)

climato_5 <- climato %>% 
  filter(cluster==5) %>% 
  select(climato)
climato_5 <- as.numeric(climato_5$climato)
climato_5 <- rep(climato_5,365)

svdMatrix_calib <- svdMatrix %>% 
  mutate(climato_cluster_1 = rep(climato_1,365),
         climato_cluster_2 = rep(climato_2,365),
         climato_cluster_3 = rep(climato_3,365),
         climato_cluster_4 = rep(climato_4,365),
         climato_cluster_5 = rep(climato_5,365))


  

weightings_clusters<-read_csv("data/weightings_clusters_all.csv") %>% select(floatID, cluster)
weightings_clusters <- weightings_clusters[which(weightings_clusters$floatID %in% colnames(svdMatrix)),]


for (i in 4:ncol(inputMatrix)){
  float <- colnames(inputMatrix)[i]
  cluster <- clusters %>% filter(floatID == float) %>% select(cluster) %>% as.numeric()
  if(cluster==1){
    inputMatrix[,i] <- inputMatrix[,i] - climato_1
  }
  if(cluster==2){
    inputMatrix[,i] <- inputMatrix[,i] - climato_2
  }
  if(cluster==3){
    inputMatrix[,i] <- inputMatrix[,i] - climato_3
  }
  if(cluster==4){
    inputMatrix[,i] <- inputMatrix[,i] - climato_4
  }
  if(cluster==5){
    inputMatrix[,i] <- inputMatrix[,i] - climato_5
  }
}

write_csv(inputMatrix,"data/svdMatrix_timeseries_DOXY_climato.csv")

  