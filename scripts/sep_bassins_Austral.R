plotData <- read_csv("data/Chla_BBP_data/plotData.csv") %>% 
  filter(cluster==1 |cluster==2 |cluster==3)

metaData <- read.csv("data/timeseries_metadata_MLD_ISO.csv") %>%
  group_by(floatID) %>%
  dplyr::summarise(floatID=unique(floatID), latitude = mean(latitude), longitude = mean(longitude)) %>% 
  ungroup()

weightings_clusters<-read_csv("data/weightings_clusters_all.csv")%>% 
  filter(cluster==1 |cluster==2 |cluster==3)

float_positions <- left_join(weightings_clusters,metaData, by = c("floatID")) %>%
  select(floatID, cluster, latitude, longitude) 

float_positions <- float_positions %>%
  dplyr::mutate(bassin = ifelse(between(longitude,-180,-65)==T,"Pacifique",ifelse(between(longitude,120,180)==T, "Pacifique", ifelse(between(longitude,20,120)==T,"Indien","Atlantique"))))
                      
#######################################

split_svdMatrix_Austral <- function(svdMatrix,float_positions,outputname){
  
  svdMatrix<- read_csv("data/svdMatrix_timeseries_env.csv")
  
  floats<-float_positions$floatID
  
  floats_Ind <- float_positions %>% filter(bassin=="Indien" & floatID %in% colnames(svdMatrix)) 
  floats_Ind <- floats_Ind$floatID
  floats_Ind <-floats_Ind #[-48]
  floats_Pac <- float_positions %>% filter(bassin=="Pacifique"& floatID %in% colnames(svdMatrix)) %>% select(floatID) 
  floats_Pac <- floats_Pac$floatID
  floats_Atl <- float_positions %>% filter(bassin=="Atlantique"& floatID %in% colnames(svdMatrix)) %>% select(floatID)
  floats_Atl <- floats_Atl$floatID
  
  floats <- floats[which(floats %in% colnames(svdMatrix))]#[-182] #5906032_2
  svdMatrix_temp <- svdMatrix %>% select(jDay,PRES,PARAM,floats)
  
  svdMatrix_Pac <- svdMatrix %>% select(jDay,PRES,PARAM,all_of(floats_Pac))
  svdMatrix_Ind <- svdMatrix %>% select(jDay,PRES,PARAM,all_of(floats_Ind))
  svdMatrix_Atl <- svdMatrix %>% select(jDay,PRES,PARAM,all_of(floats_Atl))

  # svdMatrix_final<-svdMatrix[, apply(svdMatrix, 2, function(x) !all(is.na(x)))]
  #svdMatrix_final<-cbind(svdMatrix[,1:2:3],svdMatrix_final)
  write_csv(svdMatrix,outputname)
}

split_svdMatrix(svdMatrix<- read_csv("data/Chla_BBP_data/svdMatrix_timeseries_all_variables.csv"),
                weightings_clusters<-read_csv("data/weightings_clusters_HL.csv"),
                outputname <- "data/svdMatrix_timeseries_DOXY_HL.csv")