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


###########################################################################
#réasigner les flotteurs dans clusters 1 et 2
metaData <- read.csv("data/timeseries_metadata.csv")
weightings_clusters<-read.csv("figures/vecteurs_CHL_BBP/weightings_clusters_TOT.csv")

reassign <- function(metaData,weightings_clusters,cluster_BL){
  
  library(rworldmap)
  worldMap<-getMap(resolution="coarse")
  
  # moyenner latitude et longitude de chaque TS
    lat<-vector()
    lon<-vector()
    for (i in 1:length(weightings_clusters$floatID)){
      timeserie <- which(metaData$floatID==weightings_clusters$floatID[i])
      lat[i] <- mean(metaData$latitude[timeserie])
      lon[i] <- mean(metaData$longitude[timeserie])
    }
    
    # mettre en forme la matrice
    metaData = t(matrix(data = c(weightings_clusters$floatID, lon, lat), nrow = 3, byrow = TRUE))
    colnames(metaData) <- c("floatID","longitude","latitude")
    
    # weightings.map.data<-weightings_clusters[,c(1,12)] %>%
    #   merge(metaData,by="floatID",all.x=TRUE)
    
    weightings.map.data<-weightings_clusters[,c(1,12,13)] %>%
      merge(metaData,by="floatID",all.x=TRUE)
    
    weightings.map.data$latitude<-as.numeric(weightings.map.data$latitude)
    weightings.map.data$longitude<-as.numeric(weightings.map.data$longitude)
    weightings.map.data$cluster<-as.character(weightings.map.data$cluster)
    
    # si longitude < |25| alors réassigné aux basses latitudes
    for (i in 1:length(weightings.map.data$cluster)){
      if (weightings.map.data$latitude[i] <= 25 & weightings.map.data$latitude[i] >= -25){
        weightings.map.data$cluster[i] = cluster_BL
      }
      #flotteurs Med dans les hautes lat réassignés aux basses
      if (weightings.map.data$cluster[i]!=cluster_BL & weightings.map.data$latitude[i]<=50 & weightings.map.data$latitude[i]>=25 & weightings.map.data$longitude[i]>=0 & weightings.map.data$longitude[i]<=50){
        weightings.map.data$cluster[i] = cluster_BL
      }
    }
    
    # map des TS en fonction de leur cluster
    clusterMap <- ggplot()+
      geom_polygon(data = worldMap, aes(x = long, y = lat, group=group),
                   fill = "grey", col = "black") + 
      geom_point(data = weightings.map.data,
                 aes(longitude, latitude, fill = cluster), col = "black", shape = 21, alpha = 0.8, size = 3)+
      # geom_hline(aes(yintercept = 40),linetype = 2) + 
      # geom_hline(aes(yintercept = -40),linetype =2) + 
      guides(fill=guide_legend(byrow = TRUE)) +
      theme_bw() + 
      theme(#axis.title.x = element_blank(),
            #axis.title.y = element_blank(),
            legend.title = element_blank(),
            #axis.text = element_blank(),
            #axis.ticks = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(face="bold"),
            plot.title = element_text(hjust = 0.5),
            legend.box.spacing = unit(0.1, "pt"))
    
    ggsave(plot=clusterMap,"figures/clusterAnalysis/clusterMap.png",
           width = 8, height = 5, dpi = 300)
    
    weightings_clusters$cluster<-weightings.map.data$cluster
    write_csv(weightings_clusters,"data/weightings_clusters_TOT.csv")
}



########################################################################################
#séparer matrices pour HL et BL avec cluster 1 et 2
weightings_clusters<-read_csv("data/weightings_clusters_TOT.csv")
svdWeightings<-read.csv("data/weightings_TOT.csv")

split_HL_BL <- function(weightings_clusters,svdWeightings,cluster_BL,cluster_HL){
    HL=vector()
    BL=vector()
    
    for (i in 1:length(weightings_clusters$floatID)){
      if(weightings_clusters$cluster[i]==cluster_HL){
        HL[i]<-weightings_clusters$floatID[i]
      }
      if(weightings_clusters$cluster[i]==cluster_BL){
        BL[i]<-weightings_clusters$floatID[i]
      }
    }
    HL<-HL[!is.na(HL)]
    BL<-BL[!is.na(BL)]
    
    #séparer les poids en HL et BL
    weightings_HL <- svdWeightings[-which(svdWeightings$floatID %in% BL),]
    weightings_BL <- svdWeightings[-which(svdWeightings$floatID %in% HL),]
    write_csv(weightings_HL,"data/weightings_HL.csv")
    write_csv(weightings_BL,"data/weightings_BL.csv")
}

#################################################################################
weightings_clusters<-read_csv("data/weightings_HL.csv")

split_HL <- function(weightings_clusters_HL,metaData){
  
    lat<-vector()
    lon<-vector()
    for (i in 1:length(weightings_clusters$floatID)){
      timeserie <- which(metaData$floatID==weightings_clusters$floatID[i])
      lat[i] <- mean(metaData$latitude[timeserie])
      lon[i] <- mean(metaData$longitude[timeserie])
    }
    metaData = as.data.frame(t(matrix(data = c(weightings_clusters$floatID, lon, lat), nrow = 3, byrow = TRUE)))
    colnames(metaData) <- c("floatID","longitude","latitude")
  
    floats_HLN <- metaData$floatID[which(metaData$latitude>0)]
    floats_Austral <- metaData$floatID[which(metaData$latitude<0)]
    weightings_clusters_HL_N <- weightings_clusters[which(weightings_clusters$floatID %in% floats_HLN),] 
    weightings_clusters_Austral <- weightings_clusters[which(weightings_clusters$floatID %in% floats_Austral),]
    
    write_csv(weightings_clusters_Austral,"data/weightings_Austral.csv")
    write_csv(weightings_clusters_HL_N,"data/weightings_HL_N.csv")
}

##################################################################################"
#trouver flotteurs precis
metaData <- read.csv("data/timeseries_metadata_C.csv")
weightings_clusters<-read_csv("data/all/weightings_clusters_all.csv")

found_float <- function(metaData,weightings_clusters){
  
  lat<-vector()
  lon<-vector()
  for (i in 1:length(weightings_clusters$floatID)){
    timeserie <- which(metaData$floatID==weightings_clusters$floatID[i])
    lat[i] <- mean(metaData$latitude[timeserie])
    lon[i] <- mean(metaData$longitude[timeserie])
  }
  metaData <- metaData %>% 
    filter(floatID %in% weightings_clusters$floatID) %>%
    group_by(floatID)%>%
    summarise(floatID=unique(floatID),lat=mean(latitude),long=mean(longitude))
  # metaData = as.data.frame(t(matrix(data = c(weightings_clusters$floatID, lon, lat), nrow = 3, byrow = TRUE)))
  # colnames(metaData) <- c("floatID","longitude","latitude")
  
  weightings.map.data<-weightings_clusters[,c(1,12)] %>%
    merge(metaData,by="floatID",all.x=TRUE)
  
  # weightings.map.data$latitude<-as.numeric(weightings.map.data$latitude)
  # weightings.map.data$longitude<-as.numeric(weightings.map.data$longitude)
  weightings.map.data$cluster<-as.character(weightings.map.data$cluster)
  
  # changer long, lat, cluster avec critères voulus
  TS <- weightings.map.data$floatID[which(
                                      # weightings.map.data$latitude<=-25 &
                                      # weightings.map.data$latitude>=-50 &
                                      weightings.map.data$long>=-50 &
                                      weightings.map.data$long<=0 &
                                      weightings.map.data$cluster=="4"
                                      )]
}

TS_NB <- read_csv("finalCluster_weightings.csv")
Med <- unique(TS_NB$floatID[which(TS_NB$basin == "Med_W" | TS_NB$basin == "Med_E" )])
which(TS_NB$floatID %in% Med & TS_NB$sixClusters == "SDCM")
TS_NB$sixClusters[which(TS_NB$floatID %in% Med)] 

##########################################################################
##################################################################################
#####################################################################################
#creer matrice BL et HL pour analyse cluster final

split_svdMatrix <- function(svdMatrix,weightings_clusters,outputname){
      floats<-weightings_clusters$floatID
      # BL <- vector()
      # for (i in 4:length(colnames(svdMatrix))){
      #   print(i)
      #   if(colnames(svdMatrix)[i] %in% floats == FALSE){
      #     BL[i]<- colnames(svdMatrix)[i]
      #   }
      # }
      # BL <- BL[!is.na(BL)]
      
      svdMatrix <- svdMatrix %>% select(jDay, PRES, PARAM, floats[-c(51,52,67,68,83,125,182,240,242,237)])
      # svdMatrix_final<-svdMatrix[, apply(svdMatrix, 2, function(x) !all(is.na(x)))]
      #svdMatrix_final<-cbind(svdMatrix[,1:2:3],svdMatrix_final)
      write_csv(svdMatrix,outputname)
}
    
split_svdMatrix(svdMatrix<- read_csv("data/svdMatrix_timeseries_MLD_ISO.csv"),
                weightings_clusters<-read_csv("data/weightings_clusters_all.csv"),
                outputname <- "data/svdMatrix_timeseries_MLD_ISO_HL.csv")
######################################################################################
######################################################################################
###############################################################################
  
# mld <- function(x, depth, ref.depths=5:10, criteria=c(0.03, 0.01), default.depth=NA, n.smooth=0, k=2) {
#   # check input
#   # ok <- check_input(x, depth)
#   # if (!ok) { return(NA) }
#   
#   # smooth the profile (if requested)
#   # x <- smooth(x, k=k, n=n.smooth)
#   
#   # compute the reference value
#   iref <- which(depth >= min(ref.depths) & depth <= max(ref.depths))
#   xref <- mean(x[iref], na.rm=TRUE)
#   if (is.na(xref)) {
#     warning("No data at reference depth(s).")
#     m <- NA
#   } else {
#     for (crit in criteria) {
#       i <- which(x > (xref + crit) & depth > max(ref.depths))
#       # NB: we want the element previous to meeting this criterion
#       if (length(i) > 0) {
#         a <- which(x[i]==min(x[i]))
#         i <- which(x==x[i][a])
#         break
#       }
#     }
#     
#     # extract the corresponding depth
#     m <- get_depth(i, depth)
#     
#     # replace by the default value when no criterion is met
#     if (is.na(m)) {
#       m <- default.depth
#     }
#   }
#   return(m)
# }
# 
# MLD <- function(files){
#     for(i in 1:length(files)){
#       print(paste("File",i,"of",length(files)))
#       import = read_csv(files[i])
#       currentFloat = import$float[1]
#       import$MLD = NA
#       for(j in 1:max(import$cycleNumber)){
#         if(j %in% import$cycleNumber){
#           import$MLD[which(import$cycleNumber==j)] <- mld(x = import$density[which(import$cycleNumber==j)],
#                                                           depth = import$PRES[which(import$cycleNumber==j)],
#                                                           ref.depths = 10,
#                                                           criteria = 0.03,
#                                                           default.depth=NA)
#         }
#       }
#       write_csv(import,
#                 paste("data/finalCSV/",currentFloat,".csv",sep=""))
# }  
# }


######################################################################################

# aou <- function(Sal,Temp,O2){
#   #calculate O2sat in mL/kg
#   T1 = (Temp + 273.15) / 100
#   
#   OSAT = -177.7888 + 255.5907 / T1 + 146.4813 * log(T1) - 22.2040 * T1
#   OSAT = OSAT + Sal * (-0.037362 + T1 * (0.016504 - 0.0020564 * T1))
#   OSAT = exp(OSAT)
#   
#   #convert from mL/kg to uM/kg
#   OSAT = OSAT * 1000 / 22.392
#   
#   #calculate AOU
#   AOU = OSAT - O2
#   return(AOU)
# }
# 
# 
# for(i in 1:length(files)){
#   print(paste("File",i,"of",length(files)))
#   adjImport = read_csv(files[i])
#   currentFloat = adjImport$float[1]
#   if("DOXY"%in%colnames(adjImport)){
#   adjImport$Chl_BBP = NULL
#   adjImport$AOU = oxySat(adjImport$DOXY,adjImport$TEMP,adjImport$PSAL)}
#   write_csv(adjImport,
#             paste("data/finalCSV/",currentFloat,".csv",sep=""))
#   
# }
#   

################################################################################

MLD_ISO_ZP_SIG <- function(plotData){
  MLD <- filter(plotData, PARAM == "MLD_003")
  ISO_1 <- filter(plotData, PARAM == "ISO_1")
  ZP <- filter(plotData, PARAM == "ISO_1")
  SIG <- filter(plotData, PARAM == "SIG")
  ZP$value <- pmax(ISO_1$value, MLD$value)

  write_csv(MLD,"data/plotData_MLD.csv")
  write_csv(ISO_1,"data/plotData_ISO_1.csv")
  write_csv(ZP,"data/plotData_ZP.csv")
  write_csv(SIG,"data/plotData_SIG.csv")
  
}

#######################################################################################


metaData <- read_csv("data/timeseries_metadata.csv")
#weightings_clusters_BL<-read_csv("data/weightings_clusters_BL.csv")
weightings_clusters_Austral<-read_csv("data/weightings_clusters_Austral.csv")
weightings_clusters_HL_N<-read_csv("data/weightings_clusters_HL_N.csv")
weightings_clusters_all<-read_csv("data/all/weightings_clusters_all.csv")
weightings_clusters_Austral <- weightings_clusters_all%>%filter(cluster=="1"|cluster=="2"|cluster=="3")

global_map <- function(metaData,
                       weightings_clusters){

metaData$cycleNumber <- NULL

#choose only rows with used floats
metaData <- metaData[which(metaData$floatID %in% weightings_clusters$floatID),]


#create the variables and merge them
lat<-vector()
lon<-vector()
basin <- vector()
for (i in 1:length(unique(metaData$floatID))){
  range <- which(metaData$floatID==unique(metaData$floatID)[i])
  lat[i] <- mean(metaData$latitude[range])
  lon[i] <- mean(metaData$longitude[range])
  basin[i] <- metaData$basin[range][1]
}
floatID <- unique(metaData$floatID)
metaData <- as.data.frame(cbind(floatID,basin,lat,lon))

#put the data in the order of the increasing floats
metaData <- metaData %>% arrange(floatID)

#creates new colums for the clusters and subclusters
metaData <- metaData %>%
  left_join(weightings_clusters %>% select(floatID,cluster), by = c("floatID"))
# metaData$clusters[which(metaData$floats %in% weightings_clusters$floatID)] <- "Austral"
# metaData$clusters[which(metaData$floats %in% weightings_clusters$floatID)] <- "Hautes lat Nord"
#metaData$clusters[which(metaData$floats %in% weightings_clusters_BL$floatID)] <- "Basses lat"

metaData$subclusters <- NA
metaData$subclusters[which(metaData$clusters=="Austral")] <- weightings_clusters_Austral$cluster
# metaData$subclusters[which(metaData$clusters=="Austral" & metaData$subclusters == 1)] <- 6
# metaData$subclusters[which(metaData$clusters=="Austral" & metaData$subclusters == 2)] <- 1
# metaData$subclusters[which(metaData$clusters=="Austral" & metaData$subclusters == 6)] <- 2

metaData$subclusters[which(metaData$clusters=="Hautes lat Nord")] <- weightings_clusters_HL_N$cluster
# metaData$subclusters[which(metaData$clusters=="Hautes lat Nord" & metaData$subclusters == 1)] = 4
# metaData$subclusters[which(metaData$clusters=="Hautes lat Nord" & metaData$subclusters == 2)] = 5

#metaData$subclusters[which(metaData$clusters=="Basses lat")] <- weightings_clusters_BL$cluster


metaData$lat<-as.numeric(metaData$lat)
metaData$lon<-as.numeric(metaData$lon)
# metaData$subclusters<-as.character(metaData$subclusters)

#plot the data
library(rworldmap)
library(factoextra)
library(FactoMineR)
# worldMap<-getMap(resolution="coarse")

worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]


clusterMap <- ggplot()+
  geom_polygon(data = world.df, aes(x = long, y = lat, group=group),
               fill = "grey", col = "black") + 
  geom_point(data = metaData,
             aes(lon, lat, fill = factor(cluster)), col = "black", shape = 21, alpha = 0.8, size = 3)+
  # geom_hline(aes(yintercept = 40),linetype = 2) + 
  # geom_hline(aes(yintercept = -40),linetype =2) + 
  guides(fill=guide_legend(byrow = TRUE)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(face="bold"),
        plot.title = element_text(hjust = 0.5),
        legend.box.spacing = unit(0.1, "pt")) +
  scale_fill_manual(values = c("#CC0000","#00CC00","#FFCC00","#33CCFF", "#003366"))
                    # labels = c("Austral HCB (1)","Austral MCB (2)","STZ (3)","PN (4)","AN (5)"))


ggsave(plot=clusterMap,
       paste("figures/clusterAnalysis/clusterMap.png"),
       width = 8, height = 5, dpi = 300)

rep_lat <- ggplot(data = metaData, aes(x = lat, y = subclusters, fill = subclusters)) +
  geom_boxplot() +
  theme_bw() +
  geom_vline(xintercept = -50, linetype = "dashed", color = "grey") + #SIF
  geom_vline(xintercept = -60, linetype = "dashed", color = "grey") + #PF
  geom_vline(xintercept = -40, linetype = "dashed", color = "grey") + #STF
  annotation_custom(grobTree(textGrob(expression(SIZ~~~~PAZ~~~~SAZ~~~~STZ), x=0.01,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))) +
  scale_fill_manual(name = "clusters",
                    values = c("#CC0000","#FF6666","#FFCC00","#33CCFF", "#0066CC"),
                    labels = c("Austral_HCB","Austral_MCB","Austral_LCB","PN","AN")) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  labs(x="Latitude")

ggsave(plot=rep_lat,
       paste("figures/repartition_lat.png"),
       width = 6.7, height = 3, dpi = 300)

}

global_map()

#################################################################################
#################################################################################
  
isopycnals <- function(SIG,ZP,clusters){

  plotData <- SIG
  SIGMA <- as.data.frame(matrix(NA,clusters,9))
  colnames(SIGMA)=c("cluster","ZP","Z_1000","delta","SIG1","SIG2","SIG3","SIG4","SIG5")
  
  SIG1 <- plotData
  SIG2 <- plotData
  SIG3 <- plotData
  SIG4 <- plotData
  SIG5 <- plotData

  for(i in 1:clusters){
      
      range <- which(plotData$cluster==i)
      
      ref_depth <- as.integer(max(ZP$value[range]))
      SIG_ZP <- plotData$value[range][which(plotData$PRES==ref_depth)][1]
      SIG_1000 <- mean(plotData$value[range][which(plotData$PRES==1000)],na.rm = T)
      delta <- (SIG_1000 - SIG_ZP)/6
      SIGMA$cluster[i] <- i
      SIGMA$ZP[i] <- SIG_ZP
      SIGMA$Z_1000[i] <- SIG_1000
      SIGMA$delta[i] <- delta
      SIGMA$SIG1[i] <- SIG_ZP+delta
      SIGMA$SIG2[i] <- SIG_ZP+2*delta
      SIGMA$SIG3[i] <- SIG_ZP+3*delta
      SIGMA$SIG4[i] <- SIG_ZP+4*delta
      SIGMA$SIG5[i] <- SIG_ZP+5*delta
      
      for(j in 1:length(unique(plotData$jDay))){
        #print(j)
        range_1 <- which(plotData$jDay[range]==unique(plotData$jDay)[j])
        
        SIG1$value[range][range_1] <- plotData$PRES[which(plotData$value[range][range_1] >= (SIG_ZP+delta) & plotData$value[range][range_1] < (SIG_ZP+delta+0.05))[1]]
        SIG2$value[range][range_1] <- plotData$PRES[which(plotData$value[range][range_1] >= (SIG_ZP+2*delta) & plotData$value[range][range_1] < (SIG_ZP+2*delta+0.05))[1]]
        SIG3$value[range][range_1] <- plotData$PRES[which(plotData$value[range][range_1] >= (SIG_ZP+3*delta) & plotData$value[range][range_1] < (SIG_ZP+3*delta+0.05))[1]]
        SIG4$value[range][range_1] <- plotData$PRES[which(plotData$value[range][range_1] >= (SIG_ZP+4*delta) & plotData$value[range][range_1] < (SIG_ZP+4*delta+0.05))[1]]
        SIG5$value[range][range_1] <- plotData$PRES[which(plotData$value[range][range_1] >= (SIG_ZP+5*delta) & plotData$value[range][range_1] < (SIG_ZP+5*delta+0.05))[1]]
        
      }
  }

    write_csv(SIG1,"data/plotData_SIG1.csv") 
    write_csv(SIG2,"data/plotData_SIG2.csv") 
    write_csv(SIG3,"data/plotData_SIG3.csv") 
    write_csv(SIG4,"data/plotData_SIG4.csv") 
    write_csv(SIG5,"data/plotData_SIG5.csv") 
    write_csv(SIGMA,"data/recap_isopygnes.csv")
  
 
}
##############################################################################

Depths <- function(ZP,clusters){
  
  plotData <- ZP
  plotData$PARAM <- "DEP"
  DEPTH <- as.data.frame(matrix(NA,clusters,11))
  colnames(DEPTH)=c("cluster","ZP","Z_1000","delta","DEP1","DEP2","DEP3","DEP4","DEP5","DEP6","DEP7")
  
  DEP1 <- plotData
  DEP2 <- plotData
  DEP3 <- plotData
  DEP4 <- plotData
  DEP5 <- plotData
  DEP6 <- plotData
  DEP7 <- plotData
  
  for(i in 1:clusters){
    
    for (j in 1:length(unique(import$jDay))){
    range <- which(plotData$cluster==i)
    
    ref_depth <- round(max(ZP$value[range]))
    delta <- (1000 - ref_depth)/6
    DEPTH$cluster[i] <- i
    DEPTH$ZP[i] <- ref_depth
    DEPTH$Z_1000[i] <- 1000
    DEPTH$delta[i] <- round(delta)
    DEPTH$DEP1[i] <- round(ref_depth+delta)
    DEPTH$DEP2[i] <- round(ref_depth+2*delta)
    DEPTH$DEP3[i] <- round(ref_depth+3*delta)
    DEPTH$DEP4[i] <- round(ref_depth+4*delta)
    DEPTH$DEP5[i] <- round(ref_depth+5*delta)
  
    DEP1$value[range] <- round(ref_depth+delta)
    DEP2$value[range] <- round(ref_depth+2*delta)
    DEP3$value[range] <- round(ref_depth+3*delta)
    DEP4$value[range] <- round(ref_depth+4*delta)
    DEP5$value[range] <- round(ref_depth+5*delta)
      
  }
  
  write_csv(DEP1,"data/plotData_DEP1.csv") 
  write_csv(DEP2,"data/plotData_DEP2.csv") 
  write_csv(DEP3,"data/plotData_DEP3.csv") 
  write_csv(DEP4,"data/plotData_DEP4.csv") 
  write_csv(DEP5,"data/plotData_DEP5.csv") 
  write_csv(DEPTH,"data/recap_isodepths.csv")
  
  
  }
  }
###############################################################################

