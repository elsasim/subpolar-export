
data <- read_csv("data/DOXY/svdMatrix_timeseries_DOXY.csv") %>% filter(PARAM=="AOU")
# data <- data %>% filter(PARAM=="POC_Koestner")
clusters <- read_csv("data/weightings_clusters_all.csv")

plotPrep_Elsa<-function(data,clusters,scale=FALSE){
  # output_moyenne<-data[,1:3]
  output_sd<-data[,1:3]
  for(i in 1:length(unique(clusters$cluster))){
    print(i)
    floats <- clusters %>% 
      filter(cluster==i)
    floats <- as.vector(floats$floatID)
    clusterOutput <- seq(1,nrow(output_sd))
    for(j in 1:length(unique(floats))){
      print(j)
      currentFloat = unique(floats)[j] #flotteur j du cluster
      timeSeries<-data[,grep(currentFloat, colnames(data))]
      colnames(timeSeries) <- as.character(j) #TS du flotteur
      clusterOutput <- bind_cols(clusterOutput, timeSeries) 
    }

    clusterOutput <- clusterOutput[,-1] 
    
    # moyenne <- apply(clusterOutput,1,mean)
    sd <- apply(clusterOutput,1,sd)
    # sd1 <- sd/moyenne
    
    # output_moyenne<-cbind(output_moyenne,moyenne)
    output_sd<-cbind(output_sd,sd)
    
    # colnames(output_moyenne)[(i+3)] = as.character(i)
    colnames(output_sd)[(i+3)] = as.character(i)

  }
  
  #output$month = ((output$jDay-min(output$jDay))*12)/(max(output$jDay)-min(output$jDay))
  
  # output_moyenne <-output_moyenne %>%
  #   gather(cluster,value,4:(ncol(output_moyenne))) %>%
  #   mutate(month = ((jDay-min(jDay))*12)/(max(jDay)-min(jDay)))
  
  output_sd <- output_sd %>% 
    gather(cluster,value,4:(ncol(output_sd))) %>% 
    mutate(month = ((jDay-min(jDay))*12)/(max(jDay)-min(jDay)))
  
  if(scale==TRUE){
    output <- output %>% 
      group_by(cluster,PARAM) %>% 
      mutate(value = value/sqrt(sum(value^2)/(length(value)-1)))
  }
  # write_csv(output_moyenne,"data/all/plotData/plotData_POC_Koestner.csv")
  write_csv(output_sd,"data/all/plotData_sd_AOU.csv")
}
