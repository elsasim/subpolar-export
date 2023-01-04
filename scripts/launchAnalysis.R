#The scripts requires the following packages:
# install.packages("tidyverse")
# install.packages("lubridate")
# install.packages("zoo")
# install.packages("RColorBrewer")
# install.packages("cowplot")
# install.packages("factoextra")

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


#LOAD FUNCTIONS IN ANALYSIS_FUNCTIONS.R AND 
#ACCESSORY_FUNCTIONS.R BEFORE RUNNING.
#SET WORKING DIRECTORY TO PATH OF argoSVD directory:
#setwd("/Users/nicholasbock/Projects/argo/monsoonProject/")

################################################################################
################PREPARE INTERPOLATED TIME SERIES################################
################################################################################
#interpolateData takes a list of parameters as input (e.g., "CHL","BBP","DOXY", 
#and iterates through files in the finalCSV directory to interpolate parameter 
#values horizontally and vertically. Parameter values are joined end on end
#into a single vector for each timeseries. The function saves a matrix where
#each column consists of a given timeseries.
###Arguments:###
###params: parameters to be included in profile list; only profiles containing
###all specified parameters will be retained. Default is c("BBP700","CHLA)
###tsDepth: maximum depth to be included in timeseries.
###plot: If TRUE, function will save a plot of timeseries data for each parameter

interpolateData(params=c("BBP700","CHLA_recal","ISO_1","MLD_003"),
                inputDirectory = "data/Chla_BBP_data/files_POC/",
                tsDepth=1000,plot=F)

# "TEMP","PSAL","MLD_003","ISO_1","ZP","zeta"
################################################################################
############################PERFORM SVD ANALYSIS################################
################################################################################
#svd analysis function perofms SVD on the interpolated time series produced by the
#interpolateData function, saving SVD weightings for each timeseries as a .CSV file
###Arguments:###
###svdMatrix: svdMatrix produced by interpolateData function. Default is "data/svdMatrix.csv"
###scaleWeightings: If TRUE, SVD weightings are centered and scaled before saving as .CSV.
###Default is TRUE.
###excludedFloats: Floats to be excluded from SVD analysis. Default is NA.
###plot: If TRUE, will save a scree plot for SVD results, along with plots of the first
###six singular vectors for each parameter.

#NOTE! Despite all of the rounds of QC in the import scripts, there are inevitably
#a few bad timeseries in the interpolated data. Review the timeseries figures 
#created by the interpolateData function and specify those that should be removed
#as excludedFloats. The timeseries included here are often erroneously included,
#although there may be others as well.
metaData <- read.csv("data/timeseries_metadata.csv")

################################################################################
excludeList <- c("5906043_2","6901862_3",
                 "6901485_2","6902735_2","6902735_3",
                 "6902740_2","6902907_2","6902907_3","6903574_2")

#Note that the interpolatedData function saves the interpolated data either as
#individual profiles (svdMatrx_profiles.csv), or as timeseries vectors
#"svdMatrix_timeseries.csv." The function uses "svdMatrix_timeseries" as a default;
#make sure to adjust if working with profile data
svdMatrix<- read_csv("data/svdMatrix_timeseries_env.csv")
svdMatrix_POC<- read_csv("data/svdMatrix_timeseries_POC.csv")
write_csv(svdMatrix,"data/svdMatrix_timeseries.csv")

svdAnalysis(svdMatrix = svdMatrix,
            excludedFloats = excludeList,plot = T)

################################################################################
##########################PERFORM CLUSTER ANALYSIS##############################
################################################################################
#clusterAnalysis  performs cluster analysis on  the SVD weightings produced
#by the svdAnalysis function, returning a .csv file of the classifications for
#each of the original time series, along with the corresponding weightings and
#silhouette scores (reflecting intracluster similarity. If 
#cluster parameters are not provided as function arguments, the script calls the
#clustReport function, identifying the parameters that maximize intracluster
#similarity of SVD weightings.
###Arguments:###
###weightings: data frame of weightings produced by svdAnalysis function
###clusters: optional; number of clusters to form in cluster analysis
###vectors: optional; number of SVD coefficients to use in cluster analysis
###method: clustering algorithm to use. See factoextra documentation for details
###distanceMeas: distance measure to use in clustering algorithm. See factoextra
###documentation for details
###plot: if TRUE produces a plot of vector weightings for each timeseries and cluster

svdWeightings<-read.csv("data/weightings_TOT.csv")
metaData <- read.csv("data/timeseries_metadata.csv")

clusterAnalysis(weightings = svdWeightings,
                metaData = metaData,
                vectors = 2,
                method = "diana",
                clusters = 2,
                #distanceMeas = "euclidean",
                distanceMeas = "euclidean",
                plot=TRUE,
                outputname = "weightings_clusters_TOT.csv")
#### RENAME WEIGHTING_CLUSTERS_TOT

################################################################################
############################REASSIGN############################################
################################################################################

#reasign HL and LL floats
metaData <- read.csv("data/timeseries_metadata.csv")
weightings_clusters<-read.csv("data/weightings_clusters_TOT.csv")

reassign(metaData,weightings_clusters,cluster_BL=1)

################################################################################
###########################obtain weightings for HL and LL######################
################################################################################


weightings_clusters<-read_csv("figures/vecteurs_CHL_BBP/weightings_clusters_TOT_new.csv")
svdWeightings<-read_csv("data/weightings_TOT.csv")

split_HL_BL (weightings_clusters,svdWeightings,cluster_BL=1,cluster_HL=2)

################################################################################
####################split the HL in 2 (HL_N and Austral)########################
################################################################################
weightings_clusters<-read_csv("data/weightings_clusters_all.csv")
metaData <- read_csv("data/timeseries_metadata.csv")

split_HL (weightings_clusters,metaData)

################################################################################
############## rerun the cluster analysis with BL weightings####################
################################################################################


svdWeightings<-read.csv("data/weightings_HL.csv")
metaData <- read.csv("data/timeseries_metadata.csv")

clusterAnalysis(weightings = svdWeightings,
                metaData = metaData,
                vectors = 2,
                method = "diana",
                clusters = 2,
                distanceMeas = "euclidean",
                #distanceMeas = "manhattan",
                plot=TRUE,
                outputname = "weightings_clusters_HL.csv")


################################################################################
############## rerun the cluster analysis with HL_N weightings##################
################################################################################

split_HL(weightings_clusters_HL,metaData)

svdWeightings<-read.csv("data/weightings_HL_N.csv")
metaData <- read.csv("data/timeseries_metadata.csv")

clusterAnalysis(weightings = svdWeightings,
                metaData = metaData,
                vectors = 3,
                method = "diana",
                clusters = 2,
                distanceMeas = "euclidean",
                #distanceMeas = "manhattan",
                plot=TRUE,
                outputname = "weightings_clusters_HL_N.csv")


################################################################################
############## rerun the cluster analysis with AUSTRAL weightings###############
################################################################################

svdWeightings<-read.csv("data/weightings_Austral.csv")
metaData <- read.csv("data/timeseries_metadata.csv")

clusterAnalysis(weightings = svdWeightings,
                metaData = metaData,
                vectors =2,
                method = "diana",
                clusters = 3,
                distanceMeas = "euclidean",
                # distanceMeas = "manhattan",
                plot=TRUE,
                outputname = "weightings_clusters_Austral.csv")

################################################################################
############## split matrix ####################################################
################################################################################

# Austral
split_svdMatrix(svdMatrix <- read_csv("data/svdMatrix_timeseries_Temp.csv"),
                weightings_clusters<-read_csv("data/weightings_clusters_Austral.csv"),
                outputname <- "data/svdMatrix_timeseries_Austral_POC.csv")

# HL_N
split_svdMatrix(svdMatrix <- read_csv("data/svdMatrix_timeseries_CHL_BBP.csv"),
                weightings_clusters<-read_csv("data/weightings_clusters_HL_N.csv"),
                outputname <- "data/svdMatrix_timeseries_HL_N.csv")

# BL
split_svdMatrix(svdMatrix <- read_csv("data/old/svdMatrix_timeseries_CHL_BBP.csv"),
                weightings_clusters<-read_csv("data/weightings_clusters_HL.csv"),
                outputname <- "data/svdMatrix_timeseries_HL.csv")

## All
split_svdMatrix(svdMatrix <- svdMatrix,
                weightings_clusters<-read_csv("data/weightings_clusters_all.csv"),
                outputname <- "data/svdMatrix_timeseries_all_env.csv")

################################################################################
##########################PREPARE AVG. TIMESERIES###############################
################################################################################
#plotPrep function takes the original svdMatrix and the corresponding cluster 
#analysis results as an input and averages the timeseries belonging to each cluster
#for each parameter
###Arguments###
###inputMatrix: the svdMatrix_timeseries.csv file generated by the interpolateData
###function
###clusters: the weightings_clusters.csv file produced by the clusterAnalysis function
###scale: a boolean specifying whether to scale timeseries data for each cluster 
###parameter to standard deviation of 1. If there are dramatic differences in parameter
###values between clusters, setting scale to TRUE will draw out trends more clearly
###in regions characterized by low parameter values


inputMatrix <- read_csv("data/svdMatrix_timeseries_env.csv")
inputMatrix <- svdMatrix_Atl
clusters<-read_csv("data/weightings_clusters_all.csv")

# si inputMatrix et clusters pas la même taille : remet de la même taille
which(colnames(inputMatrix) %in% clusters$floatID==F) #if result is 1,2,3 it's ok
clusters <- clusters[which(clusters$floatID %in% colnames(inputMatrix)),]

inputMatrix_temp <- inputMatrix[,which(colnames(inputMatrix) %in% clusters$floatID)]
inputMatrix <- bind_cols(inputMatrix[,(1:3)], inputMatrix_temp)
rm(inputMatrix_temp)
# créer fichier plotData (mise en nforme des données pour ploter la TS moyennée)
plotData<-plotPrep(inputMatrix,clusters %>% arrange(cluster))

# selectionner uniquement les parametres voulus (facultatif)
# plotData_POC_mask <- filter(plotData, PARAM=="POC_masked")
# write_csv(plotData_POC_mask,"data/plotData_POC_mask.csv")
# 
# plotData_POC_mask <- filter(plotData, PARAM=="POC")
# write_csv(plotData_POC_mask,"data/plotData_POC.csv")


################################################################################
###########################PLOT AVG. TIMESERIES#################################
################################################################################
#plotMaker function takes the averaged time series produced by the plotPrep
#function and creates plots as .ggplot objects, saving figures as .png files.
###Arguments###
###plotData: averaged plot data produced by plotPrep function
###ncol: number of columns in plot
###legend: if FALSE removes legend
###xAxisText: if FALSE removes text from x axis
###yAxisText: if FALSE removes text from y axis
###maxDepth: maximum depth to include on plots
###ribbon: includes standard deviation ribbon for parameters
###scaleFactor: applies exponential transformation to breaks in color gradient 

plotData<-read_csv("data/plotData.csv")

#run if matrix with chla BBP and other parameters (not with DOXY data)
MLD_ISO_ZP_SIG(plotData)  # créé les fichiers MLD, ISO, ZP correspondant aux clusters

# charger les plotData des paramètres environnementaux
MLD <- read_csv("data/plotData_MLD.csv")
ISO_1 <- read_csv("data/plotData_ISO_1.csv")
ZP <- read_csv("data/plotData_ZP.csv")

############### ANALYSIS ###############################################################
# Dall_Olmo_2014_processing ################
library(plyr)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(smoother)
library(tidyr)
library(scales)
library(grid)

# Check that it doesn't match any non-number
numbers_only <- function(x) !grepl("\\D", x)
###############################################################################
###############################################################################
### BBP
###############################################################################

#### Get data
profile_data <- read_csv("data/plotData.csv") %>% arrange(cluster)
zp_data <- read_csv("data/all/zp_data_all.csv") %>% arrange(cluster)
mld_data <- read_csv("data/all/mld_data_all.csv") %>% arrange(cluster)
mask_data <- read_csv("data/plotData.csv") %>% arrange(cluster)%>% filter(PARAM == "Cphyto_masked")
POC_data <- read_csv("data/plotData_POC_K.csv") %>% arrange(cluster) %>% filter(PARAM == "POC_Koestner")
TEMP_data <- read_csv("data/all/TEMP_data_all.csv") %>% arrange(cluster)
POC_sd <- read_csv("data/all/plotData/plotData_sd_POC.csv") %>% arrange(cluster)
AOU_sd <- read_csv("data/all/plotData/plotData_sd_AOU.csv") %>% arrange(cluster)

plotData_DOXY <- read_csv("data/all/DOXY_data_all.csv") 
plotData_SIG <- read_csv("data/all/SIG_data_all.csv")

Dall_Olmo_2014_processing_BBP (profile_data, zp_data, mld_data, mask_data, POC_data, window = 30)
  

###############################################################################
###############################################################################
## AOU
##############################################################################
profile_data_DOXY <- read_csv("data/all/DOXY_data_all.csv")
mld_data <- read_csv("data/all/mld_data_all.csv")
zp_data <- read_csv("data/all/zp_data_all.csv")
SIG_data <- read_csv("data/all/SIG_data_all.csv")
numbers_only <- function(x) !grepl("\\D", x)



################ AOU ############
Kheireddine_2020_processing_AOU (values, depths, dates, zp, depth_bin, mld, metrics)


# get data
DEP <- read_csv("data/DEP_data_all.csv")
MLD <- read_csv("data/mld_data_all.csv")
ISO_1 <- read_csv("data/ISO_data_all.csv")
ZP <- read_csv("data/zp_data_all.csv")

ZP_timings_max <- ZP %>%
  group_by(cluster) %>%
  filter(value == max(value, na.rm = T)) %>%
  summarise(max_zp_date = unique(jDay), cluster = unique(cluster)) %>%
  ungroup()

ZP_timings_min <- ZP %>%
  filter(cluster==3) %>%
  # group_by(cluster) %>%
  filter(value == min(value, na.rm = T)) %>%
  summarise(min_zp_date = unique(jDay), cluster = unique(cluster)) %>%
  ungroup()

timings <- left_join(ZP_timings_max,ZP_timings_min, by=c("cluster"))
timings <- left_join(timings,slopes_AOU, by=c("cluster"))

plotData <- read_csv("data/sep_Austral/env/plotData_Atl.csv") 

ISO_1 <- read_csv("data/sep_Austral/env/plotData_Ind.csv") %>% filter(PARAM == "ISO_1")
MLD <- read_csv("data/sep_Austral/env/plotData_Ind.csv") %>% filter(PARAM == "MLD_003")
DEP <- read_csv("data/DEP_data_all.csv")%>% filter(cluster==1|cluster==2|cluster==3)
ZP <- read_csv("data/sep_Austral/env/plotData_ZP_Ind.csv")
# charger fichier en fonction de ce qu'on veut ploter
plotData <-read_csv("data/all/mask_data_all.csv")
plotData <-read_csv("data/DOXY/plotData_climato_DOXY.csv")
plotData <-read_csv("data/all/POC_data_all.csv")
plotData <-read_csv("data/all/profile_data_all.csv")
plotData <- read_csv("data/sep_Austral/variables/plotData_Ind.csv") %>% filter(PARAM == "BBP700" | PARAM == "CHLA_recal" | PARAM == "POC_Koestner")
plotData <- read_csv("data/DOXY/plotData.csv") %>% filter(PARAM == "DOXY")


# ploter les series temporelles
plotMaker(plotData,MLD,ISO_1,ZP, DEP, params=c("MLD","ISO_1","ZP","DEP"), outputname = "") #plot with isodepths
plotMaker(plotData,MLD,ISO_1,ZP, DEP, params=c("MLD","ISO_1","ZP"), outputname = "_without_DEP") #plot without isodepths
plotMaker(plotData,MLD,ISO_1,ZP, DEP, params=NA, outputname = "_without_DEP") #plot without isodepths

################################################################################
# fonction pour moyenner les DOXY Delayed sur les clusters definis avant
################################################################################
DOXY_mean(inputMatrix <- read_csv("data/svdMatrix_timeseries_DOXY.csv"), #svdMatrix avec Chla, BBP et DOXY(D)
          clusters_BL <- read_csv("data/weightings_clusters_BL.csv"), #cluster BL de l'analyse avec chla et bbp
          clusters_HL <- read_csv("data/weightings_clusters_HL.csv")) #cluster HL de l'analyse avec chla et bbp

# si interpolate avec que DOXY (donc svdMatrix avec que param=DOXY, alors enlever
# la partie du code qui enleve les lignes Chla et BBP

################################################################################
######## GROUPER FICHIERS HAUTES ET BASSES LAT #####################
profile_data_Austral <- read_csv("data/plotData_Austral.csv")
mask_data_Austral <- read_csv("data/plotData_mask_Austral.csv")
POC_data_Austral <- read_csv("data/plotData_POC_Austral.csv")
DOXY_data_Austral <- read_csv("data/plotData_DOXY_Austral.csv")
zp_data_Austral <- read_csv("data/plotData_ZP_Austral.csv")
mld_data_Austral <- read_csv("data/plotData_MLD_Austral.csv")
ISO_data_Austral <- read_csv("data/plotData_ISO_1_Austral.csv")
DEP_data_Austral <- read_csv("data/plotData_DEP_Austral.csv")

profile_data_HLN <- read_csv("data/plotData_HLN.csv")
mask_data_HLN <- read_csv("data/plotData_mask_HLN.csv")
POC_data_HLN <- read_csv("data/plotData_POC_HLN.csv")
DOXY_data_HLN <- read_csv("data/plotData_DOXY_HLN.csv")
zp_data_HLN <- read_csv("data/plotData_ZP_HLN.csv")
mld_data_HLN <- read_csv("data/plotData_MLD_HLN.csv")
ISO_data_HLN <- read_csv("data/plotData_ISO_1_HLN.csv")
DEP_data_HLN <- read_csv("data/plotData_DEP_HLN.csv")


DOXY_data_HLN$cluster[which(DOXY_data_HLN$cluster==2)] <- 5
DOXY_data_HLN$cluster[which(DOXY_data_HLN$cluster==1)] <- 4

mask_data_HLN$cluster[which(mask_data_HLN$cluster==2)] <- 5
mask_data_HLN$cluster[which(mask_data_HLN$cluster==1)] <- 4

POC_data_HLN$cluster[which(POC_data_HLN$cluster==2)] <- 5
POC_data_HLN$cluster[which(POC_data_HLN$cluster==1)] <- 4

profile_data_HLN$cluster[which(profile_data_HLN$cluster==2)] <- 5
profile_data_HLN$cluster[which(profile_data_HLN$cluster==1)] <- 4

zp_data_HLN$cluster[which(zp_data_HLN$cluster==2)] <- 5
zp_data_HLN$cluster[which(zp_data_HLN$cluster==1)] <- 4

mld_data_HLN$cluster[which(mld_data_HLN$cluster==2)] <- 5
mld_data_HLN$cluster[which(mld_data_HLN$cluster==1)] <- 4

ISO_data_HLN$cluster[which(ISO_data_HLN$cluster==2)] <- 5
ISO_data_HLN$cluster[which(ISO_data_HLN$cluster==1)] <- 4

DEP_data_HLN$cluster[which(DEP_data_HLN$cluster==2)] <- 5
DEP_data_HLN$cluster[which(DEP_data_HLN$cluster==1)] <- 4

weightings_clusters_HL_N$cluster[which(weightings_clusters_HL_N$cluster==2)] <- 5
weightings_clusters_HL_N$cluster[which(weightings_clusters_HL_N$cluster==1)] <- 4

DOXY_data_all <- bind_rows(DOXY_data_Austral, DOXY_data_HLN)
mask_data_all <- bind_rows(mask_data_Austral, mask_data_HLN)
POC_data_all <- bind_rows(POC_data_Austral, POC_data_HLN)
profile_data_all <- bind_rows(profile_data_Austral, profile_data_HLN)
zp_data_all <- bind_rows(zp_data_Austral, zp_data_HLN)
mld_data_all <- bind_rows(mld_data_Austral, mld_data_HLN)
ISO_data_all <- bind_rows(ISO_data_Austral, ISO_data_HLN)
DEP_data_all <- bind_rows(DEP_data_Austral, DEP_data_HLN)

write_csv(DOXY_data_all,"data/all/DOXY_data_all.csv")
write_csv(mask_data_all,"data/all/mask_data_all.csv")
write_csv(POC_data_all,"data/all/POC_data_all.csv")
write_csv(profile_data_all,"data/DOXY/plotData.csv")
write_csv(zp_data_all,"data/all/zp_data_all.csv")
write_csv(mld_data_all,"data/all/mld_data_all.csv")
write_csv(ISO_data_all,"data/all/ISO_data_all.csv")
write_csv(DEP_data_all,"data/all/DEP_data_all.csv")

#inverser 1 et 2

DOXY_data_all <- read_csv("data/all/DOXY_data_all.csv")
mask_data_all <- read_csv("data/all/mask_data_all.csv")
POC_data_all <- read_csv("data/all/POC_data_all.csv")
POC_data_all <- read_csv("data/plotData.csv")
profile_data_all <- read_csv("data/DOXY/plotData.csv")
zp_data_all <- read_csv("data/all/zp_data_all.csv")
mld_data_all <- read_csv("data/all/mld_data_all.csv")
ISO_data_all <- read_csv("data/all/ISO_data_all.csv")
DEP_data_all <- read_csv("data/all/DEP_data_all.csv")


DOXY_data_all$cluster[which(DOXY_data_all$cluster==1)] <- 6
DOXY_data_all$cluster[which(DOXY_data_all$cluster==2)] <- 1
DOXY_data_all$cluster[which(DOXY_data_all$cluster==6)] <- 2

mask_data_all$cluster[which(mask_data_all$cluster==1)] <- 6
mask_data_all$cluster[which(mask_data_all$cluster==2)] <- 1
mask_data_all$cluster[which(mask_data_all$cluster==6)] <- 2

POC_data_all$cluster[which(POC_data_all$cluster==1)] <- 6
POC_data_all$cluster[which(POC_data_all$cluster==2)] <- 1
POC_data_all$cluster[which(POC_data_all$cluster==6)] <- 2

profile_data_all$cluster[which(profile_data_all$cluster==1)] <- 6
profile_data_all$cluster[which(profile_data_all$cluster==2)] <- 1
profile_data_all$cluster[which(profile_data_all$cluster==6)] <- 2

zp_data_all$cluster[which(zp_data_all$cluster==1)] <- 6
zp_data_all$cluster[which(zp_data_all$cluster==2)] <- 1
zp_data_all$cluster[which(zp_data_all$cluster==6)] <- 2

mld_data_all$cluster[which(mld_data_all$cluster==1)] <- 6
mld_data_all$cluster[which(mld_data_all$cluster==2)] <- 1
mld_data_all$cluster[which(mld_data_all$cluster==6)] <- 2

ISO_data_all$cluster[which(ISO_data_all$cluster==1)] <- 6
ISO_data_all$cluster[which(ISO_data_all$cluster==2)] <- 1
ISO_data_all$cluster[which(ISO_data_all$cluster==6)] <- 2

DEP_data_all$cluster[which(DEP_data_all$cluster==1)] <- 6
DEP_data_all$cluster[which(DEP_data_all$cluster==2)] <- 1
DEP_data_all$cluster[which(DEP_data_all$cluster==6)] <- 2

plotData$cluster[which(plotData$cluster==1)] <- 6
plotData$cluster[which(plotData$cluster==2)] <- 1
plotData$cluster[which(plotData$cluster==6)] <- 2

weightings_clusters_Austral$cluster[which(weightings_clusters_Austral$cluster==1)] <- 6
weightings_clusters_Austral$cluster[which(weightings_clusters_Austral$cluster==2)] <- 1
weightings_clusters_Austral$cluster[which(weightings_clusters_Austral$cluster==6)] <- 2

weightings_clusters_all <- bind_rows(weightings_clusters_Austral, weightings_clusters_HL_N)
write_csv(weightings_clusters_all,"data/all/weightings_clusters_all.csv")
write_csv(mask_data_all,"data/all/mask_data_all.csv")

