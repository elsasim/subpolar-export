################################################################################
################################################################################
################################################################################
#interpolateData function takes a directory of BGC-Argo .csv files and returns a matrices of CHLA and BBP700
#time series interpolated to 1 meter and 1 day intervals for SVD analysis. Input files must contain the following columns:
#float = Float number
#date (determined from JULD + REFERENCE_DATE)
#PRES
#CHLA
#BBP700
#tsDepth: determines the maximum depth to be included for each time series
#plot: if TRUE, function will save plots of interpolated CHLA and BBP values to local directory (as specified on line
################################################################################
################################################################################
grep("2902196",list.files(inputDirectory))
interpolateData <- function(inputDirectory = "data/finalCSV/",
                            matrixType = "timeseries",
                            params = c("CHLA","BBP700","MLD_003","ISO_1"),
                            tsDepth=1000,plot = FALSE){
  require(zoo)
  require(lubridate)
  require(tidyverse)
  require(RColorBrewer)
  require(cowplot)
  require(patchwork)

  files<-list.files(inputDirectory,full.names = T,pattern=".csv")
  
  tsTemplate<-data.frame()
  for(i in 1:length(params)){
    temp<-data.frame(jDay = rep(seq(1,365,1),each = tsDepth + 1),
                     PRES = rep(seq(0,tsDepth,1),365),
                     PARAM = params[i])
    tsTemplate<-rbind(tsTemplate,temp)
  }
  
  tsOutput<-tsTemplate %>% 
    arrange(jDay,PRES,PARAM)
  
  profOutput <- data.frame(PRES = rep(seq(0,tsDepth,1),times=length(params)),
                              PARAM = rep(params,each=1001))

  metaData <- data.frame()

  for(j in 1:length(files)){
    print(paste("File ", j, " of ", length(files),sep=""))
    inputData<-read_csv(files[j],show_col_types = F)

    #Skip file if any of specified parameters are not present in float data
    if(any(params %in% colnames(inputData)==FALSE)){
      next}
    
    metaTemp <- inputData %>%
      dplyr::select(basin,floatID,cycleNumber,date,latitude,longitude) %>%
      group_by(basin,floatID,cycleNumber,date) %>%
      summarize_all(median,na.rm=TRUE)
    
    inputData<-inputData%>% 
      dplyr::select(float,floatID,jDay,cycleNumber,PRES,all_of(params)) %>% 
      arrange(floatID,jDay,PRES) 
    
  if(nrow(inputData)==0){
    next}
    
  #Iteratively subset original float data to only include measurements from an individual year
  for(k in 1:length(unique(inputData$floatID))){
    print(paste("Time series",k,"of",length(unique(inputData$floatID))))

      #Remove missing CHLA or BBP700 values, along with any profiles where all measurements < 300 m are 0
      inputSubset = subset(inputData,floatID==unique(inputData$floatID)[k]) %>% 
        gather(PARAM,VALUE,params) %>% 
        group_by(jDay,PARAM) %>% 
        filter(all(is.na(VALUE))==FALSE,
               ####added the following line to only include the first cycle
               ####on dates where there are more than 1 cycles
               cycleNumber == min(cycleNumber,na.rm=TRUE)) %>% 
        ungroup()
      
    if(matrixType == "timeseries"){
        if(length(unique(inputSubset$PARAM))<length(params)){
          print("Missing parameters")
          next}
        if (length(unique(inputSubset$jDay[which(inputSubset$PARAM == "ISO_1")]))==1 & "ISO_1" %in% params){ # if ISO_1 only for 1 day in the year skip the TS
          print("Not enough ISO data")
          next
        }
        #Create vector of Julian dates within time series
        uniqueDays <- unique(inputSubset$jDay)
        #Calculate maximum interval between profiles 
        maxInt <- max(diff(uniqueDays),na.rm=TRUE)
        
        #Skip any time series where the first profile is completed after JDAY = 20 or where the
        #last profile is completed before JDAY = 345. Also skip any time series where there is 
        #an interval of more than 20 days between any two profiles.
        if(min(uniqueDays,na.rm=TRUE) > 20 | max(uniqueDays,na.rm=TRUE) < 345 | maxInt > 20){
          print("Insufficient profile resolution")
          next
        }
        
        inputSubset.interp<-inputSubset %>%
          dplyr::select(jDay,PRES,PARAM,VALUE) %>%
          #Merge existing data to JDAY * PRES template
          merge(tsTemplate, by = c("jDay","PRES","PARAM"),all.y=TRUE) %>%
          #Vertically interpolate CHLA and BBP700 values
          dplyr::group_by(jDay,PARAM) %>%
          dplyr::mutate(VALUE = na.fill(VALUE,fill="extend")) %>%
          ungroup() %>%
          #Horizontally interpolate CHLA and BBP700 values
          dplyr::group_by(PRES,PARAM) %>%
          dplyr::mutate(VALUE = na.fill(VALUE,fill="extend")) %>% 
          dplyr::arrange(jDay,PRES,PARAM)
        
        tsOutput <- cbind(tsOutput,inputSubset.interp$VALUE)
        colnames(tsOutput)[ncol(tsOutput)] = unique(inputData$floatID)[k]
      }
      
    if(matrixType == "profiles"){
      profTemplate<-data.frame()
      profID = paste(inputSubset$float[1],unique(inputSubset$cycleNumber),sep="_")
      for(i in 1:length(params)){
        temp<-data.frame(profID = rep(profID,each = tsDepth + 1),
                         PRES = rep(seq(0,tsDepth,1),times = length(profID)),
                         PARAM = params[i])
        profTemplate<-rbind(profTemplate,temp)
        profTemplate<-profTemplate %>% 
          arrange(profID,PRES,PARAM)
      }
      
      inputSubset.interp<-inputSubset %>%
        mutate(profID = paste(float,cycleNumber,sep="_")) %>% 
        dplyr::select(profID,PRES,PARAM,VALUE) %>%
        #Merge existing data to JDAY * PRES template
        merge(profTemplate, by = c("profID","PRES","PARAM"),all.y=TRUE) %>%
        #Vertically interpolate CHLA and BBP700 values
        group_by(profID,PARAM) %>%
        mutate(VALUE = na.fill(VALUE,fill="extend")) %>%
        pivot_wider(names_from = profID,values_from=VALUE) %>% 
        arrange(PARAM,PRES)
      profOutput <- merge(profOutput,inputSubset.interp,by=c("PRES","PARAM"))
    }
      
    if(plot==TRUE & matrixType == "timeseries"){
      outputPlots <- list()
      plotSubset<-tsOutput[,c(1:3,ncol(tsOutput))]
      colnames(plotSubset)[ncol(plotSubset)] = "VALUE"
      for(l in 1:2){
        
      outputPlots[[l]] <- ggplot() +
        scale_x_reverse(expand=c(0,0),limits=c(200,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        coord_flip() + 
        geom_tile(data = plotSubset[plotSubset$PARAM==params[l],],
                  aes(x = PRES,
                      y = jDay,
                      fill = VALUE)) +
        scale_fill_distiller(palette = "Spectral") + 
        theme(legend.title = element_blank(),
              legend.text = element_text(size =8),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill = NA),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              #axis.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5,size = 10,face="bold")) +
        labs(x = "", y = " ") +
        ggtitle(params[l])
      
      if("MLD_003" %in% params){
        outputPlots[[l]] = outputPlots[[l]] + 
          geom_path(data = plotSubset[plotSubset$PARAM=="MLD_003",], aes(x = VALUE, y = jDay),col="black",linetype = 2)
      }
      
      if("ISO_1" %in% params){
        outputPlots[[l]] = outputPlots[[l]] + 
          geom_path(data = plotSubset[plotSubset$PARAM=="ISO_1",], aes(x = VALUE, y = jDay),col="grey",linetype = 2)
      }
      
      }
      rm(plotSubset)
      finalPlot <- wrap_plots(outputPlots,nrow=2)
      ggsave(plot = finalPlot,
             paste("figures/interpolatedTimeseries/O2/PN/",unique(inputData$floatID)[k],".png",sep=""), width = 4, height = 3,
             dpi = 300)
    }
  }
    #Add current float data to metadata file
    metaData <- rbind(metaData,metaTemp)
  }
  
  #Return list containing matrices of interpolated CHLA and BBP data.
  #Can be replaced with write.csv command to save matrices to a local directory
  if(matrixType=="timeseries"){
    write_csv(tsOutput,
            "data/svdMatrix_timeseries_env.csv")

    metaData <- metaData %>%
      dplyr::select(basin,floatID,cycleNumber,latitude,longitude) %>%
      summarize(latitude = median(latitude,na.rm=TRUE),
                longitude = median(longitude,na.rm=TRUE))

    write_csv(metaData,
              "data/timeseries_metadata_ISO_MLD.csv")
  }

  if(matrixType=="profiles"){
  write_csv(profOutput,
            "data/svdMatrix_profiles.csv")

    write_csv(metaData,
              "data/profile_metadata.csv")
  }

}
################################################################################
################################################################################
################################################################################
#svd analysis function perofms SVD on the interpolated time series produced by the
#interpolateData function, returning a scree plot, plots of the first six 
#singular vectors for chl and bbp, and a .csv file of the SVD weightings
################################################################################
################################################################################

svdAnalysis<-function(svdMatrix = "data/svdMatrix_timeseries.csv",
                      excludedFloats = NA,
                      scaleWeightings = TRUE,
                      plot = FALSE){
  
  # require(vroom)
  # svdMatrix <- vroom::vroom(svdMatrix) %>% 
  #   arrange(PARAM,jDay,PRES)
  svdMatrix <- svdMatrix %>%
    filter(PRES<=750)
  
  template<-svdMatrix %>% 
    select(PARAM,jDay,PRES)
  
  params = unique(svdMatrix$PARAM)
  #colnames(svdMatrix)[4:ncol(svdMatrix)]<-substring(colnames(svdMatrix[4:ncol(svdMatrix)]),2,10)
  if(all(is.na(excludedFloats))==FALSE){
    svdMatrix<-svdMatrix[,-which(colnames(svdMatrix)%in%excludedFloats)]
  }
  
  floatIDs<-colnames(svdMatrix)[4:ncol(svdMatrix)]
  #Convert svdMatrix to matrix and remove JDAY, PRES, and PARAM columns
  svdMatrix.c<-as.matrix(svdMatrix[,4:ncol(svdMatrix)])
  
  #Center and scale each column of matrix for each parameter
  for(i in unique(svdMatrix$PARAM)){
    print(i)
    svdMatrix.c[svdMatrix$PARAM==i,] = apply(svdMatrix.c[svdMatrix$PARAM==i,],2,scale,center=TRUE,scale=TRUE)
  }
  
  #Transpose matrix
  svdMatrix.ct = t(svdMatrix.c)
  rm(svdMatrix.c)
  
  #Preform svd on final matrix
  svdResults <- svd(svdMatrix.ct)
  rm(svdMatrix.ct)
  
  #Prepare scree plot (shows percentage of variance accounted for by each eigenvector,
  #based on relative magnitude of corresponding eigenvalue (stored in svdResults$d matrix)
  totalMatrix.EV <- (svdResults$d^2)/sum(svdResults$d^2) * 100
  totalMatrix.EV[1:6]
  if(plot==TRUE){
    screePlot <- ggplot() +
      geom_line(aes(x = seq(1:length(totalMatrix.EV)),
                    y = totalMatrix.EV))+
      geom_point(aes(x = seq(1:length(totalMatrix.EV)),
                     y = totalMatrix.EV),shape = 21, size = 2, fill = "red")+
      theme_bw() +
      labs(x = "Singular Vector",
           y = "% Variance") +
      theme(plot.title = element_text(hjust = 0.5))+
      ggtitle(paste("SV1 - SV6 capture",
                    round(sum(totalMatrix.EV[1:6])/sum(totalMatrix.EV)*100,2),
                    "% of total variance"))
    
    ggsave("figures/screePlot.png",
           plot=screePlot,dpi=300)
    
    for(j in 1:length(params)){
      tempPlot <- vectorPlot(svdResults,params[j],template)
      ggsave(plot=tempPlot,
             paste("figures/",params[j],"_vectors.png",sep=""))
    }
    
    #Calculate SVD coefficients as product of left singular vectors and eigenvalues
    weightings = sweep(svdResults$u, MARGIN=2, svdResults$d,"*")
    colnames(weightings) <- paste("V_",seq(1,ncol(weightings),1),sep="")
    weightings<-data.frame(cbind(floatID = floatIDs,weightings))
    weightings[2:ncol(weightings)] <- apply(weightings[2:ncol(weightings)],
                                            2,as.numeric)
    if(scaleWeightings==TRUE){
      weightings[,2:ncol(weightings)] = apply(weightings[,2:ncol(weightings)],2,scale,scale=TRUE,center=FALSE)
    }
    
    
    rownames(weightings) = weightings$floatID
    
    #Return list of SVD weightings and corresponding vectors
    write_csv(weightings[,1:11],
              "data/weightings_TOT.csv")
  }
}
################################################################################
################################################################################
################################################################################
#vectorPlot function is called by the svdAnalysis function to produce time series
#of the first six singular vectors for chl and bbp
#svdData: svdResults provided by svd function
#param: "CHL" or "BBP"; determines which singular vectors to plot
#template: a two-column dataframe containing PRES and JDAY values for plot
################################################################################
################################################################################

vectorPlot<-function(svdData,param,template){
  require(cowplot)
  plotList<-list()
  
  svdData.EV <- (svdData$d^2)/sum(svdData$d^2) * 100
  for(i in 1:6){
    perVar<-round(sum(svdData.EV[i])/sum(svdData.EV),2)*100
    vectors<-data.frame(value=svdData$v[,i])
    #Binding the vectors to chlMetaCols prevents all 600 rows of the data from being plotted
    vectors<-cbind(template,vectors)
    plotList[[i]]<-ggplot(data = subset(vectors,PARAM==param),
                          aes(x=PRES,y=jDay,fill = value)) +
      scale_x_reverse(expand=c(0,0),limits=c(300,0)) + 
      scale_y_continuous(expand=c(0,0)) +
      coord_flip() + 
      geom_tile()+
      scale_fill_distiller(palette = "Spectral") + 
      theme(legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            #legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill = NA),
            plot.title = element_text(hjust = 0.5,size = 10, face = "bold"),
            strip.text = element_text(face="bold",size=12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12))+
      labs(x = "Depth", y = "Julian Date") +
      ggtitle(paste("SV ",i,": ",perVar,"% of variability",sep=""))
  }
  
  vectorPlot_output<-plot_grid(plotList[[1]],plotList[[2]],
                               plotList[[3]],plotList[[4]],
                               plotList[[5]],plotList[[6]],
                               ncol=2)
  return(vectorPlot_output)
}

################################################################################
################################################################################
################################################################################
#The clusterAnalysis function performs cluster analysis on  the SVD weightings produced
#by the svdAnalysis function, returning a .csv file of the classifications for
#each of the original time series, along with the corresponding weightings and
#silhouette scores (reflecting intracluster similarity. If 
#cluster parameters are not provided as function arguments, the script calls the
#clustReport function, identifying the parameters that maximize intracluster
#similarity of SVD weightings.
#weightings: data frame of weightings produced by svdAnalysis function
#clusters: optional; number of clusters to form in cluster analysis
#vectors: optional; number of SVD coefficients to use in cluster analysis
#method: clustering algorithm to use. See factoextra documentation for details
#distanceMeas: distance measure to use in clustering algorithm. See factoextra
#documentation for details
#plot: if TRUE produces a plot of vector weightings for each timeseries and cluster
################################################################################
################################################################################

clusterAnalysis <- function(weightings,
                            metaData,
                            clusters=NA,
                            vectors = NA, 
                            method = "diana",
                            distanceMeas=NA,
                            plot=FALSE,
                            outputname = "weightings_clusters.csv"){
  
  library(rworldmap)
  library(factoextra)
  library(FactoMineR)
  worldMap<-getMap(resolution="coarse")
  rownames(weightings) <- weightings$floatID
  weightings<-weightings[,-1]
  
  #If clusters, number of vectors to use in clustering, or the distance measure
  #to be use are not specified, the function will run clusterTest to iteratively
  #evaluate silhouette scores for different combinations of these factors, saving
  #the output as "data/clusterReport.csv." 
  if(is.na(clusters)|is.na(distanceMeas)|is.na(vectors)){
    clusterTest <- clustReport(weightings)
    write_csv(clusterTest,
              "data/clusterReport.csv")
    clusters = clusterTest$clusters[clusterTest$avgSil==max(clusterTest$avgSil,na.rm=TRUE)]
    distanceMeas = clusterTest$dist[clusterTest$avgSil==max(clusterTest$avgSil,na.rm=TRUE)]
    vectors = clusterTest$lastVector[clusterTest$avgSil==max(clusterTest$avgSil,na.rm=TRUE)]
  }
  
  clustResult<-eclust(weightings[,c(1:vectors)],
                      FUNcluster = "diana",
                      hc_metric = distanceMeas,
                      stand = T,
                      k = clusters, nstart =100)
  silResult<-fviz_silhouette(clustResult)
  
  ggsave(plot=silResult,
         paste("figures/clusterAnalysis/silhouette_",clusters,"_clusters_",vectors,"_vectors.png",sep=""),
         width = 5, height = 5, dpi = 300)
  
  clustResult<-data.frame(floatID = rownames(clustResult$silinfo$widths),
                          cluster = clustResult$silinfo$widths$cluster,
                          silScore = clustResult$silinfo$widths$sil_width)
  
  weightings$floatID <- rownames(weightings)
  weightings<-merge(weightings,clustResult,by="floatID")
  weightings<-left_join(weightings,clustResult,by="floatID")
  write_csv(weightings,
            paste("data/",outputname,sep=""))
  
# Moyenner les lat/lon
  lat<-vector()
  lon<-vector()
  for (i in 1:length(svdWeightings$floatID)){
      timeserie <- which(metaData$floatID==svdWeightings$floatID[i])
      lat[i] <- mean(metaData$latitude[timeserie])
      lon[i] <- mean(metaData$longitude[timeserie])
  }
  metaData = as.data.frame(t(matrix(data = c(svdWeightings$floatID, lon, lat), nrow = 3, byrow = TRUE)))
  colnames(metaData) <- c("floatID","longitude","latitude")
  
  if(plot==TRUE){
    weightings.map.data<-clustResult %>%
      merge(metaData,by="floatID",all.x=TRUE)

    weightings.plot.data<-weightings %>%
      gather(vectors,value,2:7)
    
    weightings.map.data$latitude<-as.numeric(weightings.map.data$latitude)
    weightings.map.data$longitude<-as.numeric(weightings.map.data$longitude)
    weightings.map.data$cluster<-as.character(weightings.map.data$cluster)

    clusterMap <- ggplot()+
      geom_polygon(data = worldMap, aes(x = long, y = lat, group=group),
                   fill = "grey", col = "black") + 
      geom_point(data = weightings.map.data,
                 aes(longitude, latitude, fill = cluster), col = "black", shape = 21, alpha = 0.8, size = 3)+
      # geom_hline(aes(yintercept = 40),linetype = 2) + 
      # geom_hline(aes(yintercept = -40),linetype =2) + 
      guides(fill=guide_legend(byrow = TRUE)) +
      theme_bw() + 
      theme(#axis.title.x = element_text("longitude"),
            #axis.title.y = element_text("latitude"),
            legend.title = element_blank(),
            #axis.text = element_blank(),
            #axis.ticks = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(face="bold"),
            plot.title = element_text(hjust = 0.5),
            legend.box.spacing = unit(0.1, "pt"))
    
    ggsave(plot=clusterMap,
           paste("figures/clusterAnalysis/clusterMap_",clusters,"_clusters_",vectors,"_vectors.png",sep=""),
           width = 8, height = 5, dpi = 300)
    
    weightings.plot.data$vector <- factor(weightings.plot.data$vector,
                                          levels = unique(weightings.plot.data$vector),
                                          labels = seq(1:length(unique(weightings.plot.data$vector))))
    
    weightings.plot<-ggplot(weightings.plot.data) + 
      geom_path(aes(x = as.numeric(vector) , y = value,group= floatID)) +  
      facet_wrap(~cluster) + 
      theme_bw() + 
      theme(strip.background = element_blank(),
            strip.text = element_text(face="bold")) + 
      labs(x = "Vector", y = "value")
    
    ggsave(plot = weightings.plot,
           paste("figures/clusterAnalysis/weightings_",clusters,"_clusters_",vectors,"_vectors.png",sep=""),
           width = 12, height = 5, dpi = 300)
  }
}

################################################################################
################################################################################
################################################################################
#séparer matrices pour HL et BL avec cluster 1 et 2
# HL=vector()
# BL=vector()
# weightings_clusters<-read_csv("data/weightings_clusters_TOT.csv")
# svdWeightings<-read_csv("data/weightings_TOT.csv")
# for (i in 1:length(weightings_clusters$floatID)){
#   if(weightings_clusters$cluster[i]==1){
#     HL[i]<-weightings_clusters$floatID[i]
#   }
#   if(weightings_clusters$cluster[i]==2){
#     BL[i]<-weightings_clusters$floatID[i]
#   }
# }
# HL<-HL[!is.na(HL)]
# BL<-BL[!is.na(BL)]
# 
# 
# #séparer les poids en HL et BL
# weightings_HL <- svdWeightings[-which(svdWeightings$floatID %in% BL),]
# weightings_BL <- svdWeightings[-which(svdWeightings$floatID %in% HL),]
# write_csv(weightings_HL,"data/weightings_HL.csv")
# write_csv(weightings_BL,"data/weightings_BL.csv")
################################################################################
################################################################################
################################################################################
#The clustReport function is called by cluster analysis to iteratively perform 
#cluster analysis on the SVD weightings, using varying numbers of clusters, weightings,
#and distance metrics. The function returns a data frame of silhouette scores
#for different parameter settings.
################################################################################
################################################################################

clustReport<-function(weightings){
  clusterList<-list()
  distList <- list()
  metricList <- list()
  lastVectorList <- list()
  silAvgList <- list()
  silMaxList <- list()
  silMinList <- list()
  minSizeList <- list()
  #maxVector = 
  #for(vectors in 2:min(5,ncol(weightings))){
  for(vectors in 2:6){
    for(clusters in 2:6){
      for(metric in c("euclidean","manhattan")){
        print(paste(1,"to",vectors))
        print(metric)
        clustResult<-eclust(weightings[,c(1:vectors)],
                            FUNcluster = "diana",
                            hc_metric = metric,
                            stand = F,
                            k = clusters, nstart =25)
        
        silGroup <- data.frame(clustResult$silinfo$widths) %>% 
          group_by(cluster) %>% 
          dplyr::summarize(avgSil = mean(sil_width),n = n())
        
        lastVectorList <- c(lastVectorList,vectors)
        clusterList <- c(clusterList,clusters)
        distList <- c(distList,metric)
        silAvgList <- c(silAvgList,round(clustResult$silinfo$avg.width,2))
        silMaxList <- c(silMaxList, round(max(clustResult$silinfo$clus.avg.widths),2))
        silMinList <- c(silMinList, round(min(clustResult$silinfo$clus.avg.widths),2))
        minSizeList <- c(minSizeList, min(silGroup$n))
      }
    }
  }
  
  
  report<-data.frame(lastVector = unlist(lastVectorList),
                     clusters = unlist(clusterList),
                     dist = unlist(distList),
                     minSil = unlist(silMinList),
                     maxSil = unlist(silMaxList),
                     avgSil = unlist(silAvgList),
                     minSize = unlist(minSizeList)) %>% 
    filter(minSil > 0)
  
  return(report)
}

################################################################################
################################################################################
################################################################################
#The plotPrep function takes the interpolated time series produced by the
#interpolateData function, and averages the measurements based on the .csv
#file provided by the clusterAnalysis function.
#data: interpolated time series produced by the interpolateData function
#clusters: data frame containing floatID and cluster column
#scale: if true, will scale data in each cluster by the RMS
################################################################################
################################################################################

plotPrep <- function(data,clusters,scale=FALSE){
  output <- data[,1:3]
  for(i in 1:length(unique(clusters$cluster))){
    print(i)
    floats <- clusters$floatID[clusters$cluster==unique(clusters$cluster)[i]]
    clusterOutput <- rep(0,nrow(output))
    for(j in 1:length(unique(floats))){
      print(j)
      currentFloat = unique(floats)[j] #flotteur j du cluster
      timeSeries<-data[,grep(currentFloat, colnames(data))] #TS du flotteur
      clusterOutput = clusterOutput + timeSeries
    }
    clusterOutput = clusterOutput/length(unique(floats))
    #clusterOutput = clusterOutput/length(unique(floats)[j])
    output<-cbind(output,clusterOutput)
    colnames(output)[(i+3)] = as.character(i)
  }
  
  #output$month = ((output$jDay-min(output$jDay))*12)/(max(output$jDay)-min(output$jDay))
  
  output<-output %>% 
    gather(cluster,value,4:(ncol(output))) %>% 
    mutate(month = ((jDay-min(jDay))*12)/(max(jDay)-min(jDay)))
  
  if(scale==TRUE){
    output <- output %>% 
      group_by(cluster,PARAM) %>% 
      mutate(value = value/sqrt(sum(value^2)/(length(value)-1)))
  }
  write_csv(output,"data/plotData.csv")
}
################################################################################
################################################################################
################################################################################
#The plotMaker function takes the averaged time series produced by the plotPrep
#function and creates plots as .ggplot objects, saving figures as .png files.
#plotData: averaged plot data produced by plotPrep function
#params: environmental parameters to include in plot. Options include MLD, NCLINE, ISO_01
#ncol: number of columns in plot
#legend: if FALSE removes legend
#xAxisText: if FALSE removes text from x axis
#yAxisText: if FALSE removes text from y axis
#maxDepth: maximum depth to include on plots
#ribbon: includes standard deviation ribbon for parameters
#scaleFactor: applies exponential transformation to breaks in color gradient 
################################################################################
################################################################################

plotMaker <- function(plotData,
                      MLD,ISO_1,ZP,DEP,params=NA,
                      outputname, ncol=2,
                      legend=TRUE,xAxisText=TRUE,
                      yAxisText=TRUE,maxDepth=1000,
                      ribbon=FALSE,scaleFactor = 1){
  
  if("DEP" %in% params){
  DEP1 <- dplyr::select(DEP,value = depth_0, month = month, cluster = cluster)
  DEP2 <- dplyr::select(DEP,value = depth_100, month = month, cluster = cluster)
  DEP3 <- dplyr::select(DEP,value = depth_200, month = month, cluster = cluster)
  DEP4 <- dplyr::select(DEP,value = depth_300, month = month, cluster = cluster)
  DEP5 <- dplyr::select(DEP,value = depth_400, month = month, cluster = cluster)
  DEP6 <- dplyr::select(DEP,value = depth_500, month = month, cluster = cluster)
  DEP7 <- dplyr::select(DEP,value = depth_600, month = month, cluster = cluster)
  DEP8 <- dplyr::select(DEP,value = depth_700, month = month, cluster = cluster)
  }
  
  
  for(i in 1:length(unique(plotData$PARAM))){
    plotTemp <- plotData[plotData$PARAM == unique(plotData$PARAM)[i],]
    
    avgPlot<-ggplot() +
      scale_x_reverse(expand=c(0,0),limits=c(1000,0)) +
      scale_y_continuous(expand=c(0,0)) +
      coord_flip() +
      geom_tile(data = plotTemp,
                aes(x = PRES,
                    y = month,
                    fill = value)) +
      scale_fill_distiller(palette = "Spectral",
                          # breaks = c(-5,0,5),
                          # limits=c(-5,5),
                          # oob = scales::squish,
                          # labels = c("<-5","0",">5"),
                          # name = expression(DOXY~anomaly(µmol~kg^-1~d^1))
                           ) +
      facet_wrap(~factor(cluster),ncol=3) +
      theme(legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            plot.title = element_text(hjust = 0.5,face="bold",size=12),
            strip.text = element_text(face="bold",size=12),
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 12),
            #axis.title = element_blank(),
            #axis.ticks.x = element_blank(),
            panel.border = element_rect(fill = NA)) +
      labs(x = "Depth (m)", y = "relative Month",
           fill = unique(plotData$PARAM)[i])
    
    # avgPlot<-ggplot() +
    #   scale_x_reverse(expand=c(0,0),limits=c(1000,0)) + 
    #   scale_y_continuous(expand=c(0,0)) +
    #   coord_flip() + 
    #   geom_tile(data = plotTemp,
    #             aes(x = PRES,
    #                 y = month,
    #                 fill = value)) +
    #   scale_fill_gradient2(midpoint = 0,
    #                         low = "blue",
    #                         mid = "white",
    #                         high = "red",
    #                         space = "Lab",
    #                        breaks = c(-0.1,0,0.1),
    #                        limits=c(-0.1,0.1),
    #                        oob = scales::squish,
    #                        labels = c("<-0.1","0",">0.1"),
    #                        name = expression(O2~rates(µmol~kg^-1~d^1))
    #   ) +
    #   facet_wrap(~factor(cluster),ncol=3) +
    #   theme(legend.title = element_text(size = 10),
    #         legend.text = element_text(size = 10),
    #         panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(),
    #         strip.background = element_blank(),
    #         plot.title = element_text(hjust = 0.5,face="bold",size=12),
    #         strip.text = element_text(face="bold",size=12),
    #         axis.text = element_text(size = 10),
    #         axis.title = element_text(size = 12),
    #         #axis.title = element_blank(),
    #         #axis.ticks.x = element_blank(),
    #         panel.border = element_rect(fill = NA)) +
    #   labs(x = "Depth (m)", y = "relative Month",
    #        fill = unique(plotData$PARAM)[i])
    
    if(legend==FALSE){
      avgPlot = avgPlot +
        theme(legend.title = element_blank(),
              legend.position="none")}
    
    if(xAxisText==FALSE){
      avgPlot = avgPlot +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position="none")}
    
    
    if("ISO_1" %in% params){
      avgPlot = avgPlot + 
        #geom_path(data = plotData[[2]], aes(x = ISO_10.avg, y =month, group = cluster),col="white",linetype = 2)+
        geom_path(data = ISO_1, aes(x = value, y =month, group = cluster),col="red",linetype = 2)
      #geom_path(data = plotData[[2]], aes(x = ISO_01.avg, y =month, group = cluster),col="black",linetype = 2)
    }
    
    if("IPAR20" %in% params){
      avgPlot = avgPlot + 
        #geom_path(data = plotData[[2]], aes(x = ISO_10.avg, y =month, group = cluster),col="white",linetype = 2)+
        #geom_path(data = plotData[[2]], aes(x = ISO_1.avg, y =month, group = cluster),col="grey",linetype = 2)
        geom_path(data = plotData[[2]], aes(x = IPAR20.avg, y =month, group = cluster),col="white",linetype = 2)
    }
    
    if("MLD" %in% params){
      avgPlot = avgPlot + 
        geom_path(data = MLD, aes(x = value, y =month, group = cluster),col="black",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }
    
    if("ZP" %in% params){
      avgPlot = avgPlot + 
        geom_path(data = ZP, aes(x = value, y =month, group = cluster),col="green",linetype = 4)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }
    
    # if("DEP" %in% params){
    #   avgPlot = avgPlot +
    #     geom_path(data = DEP, aes(x = value, y =jDay, group = cluster),col="grey",linetype = 2)
    #   if(ribbon == TRUE){
    #     avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
    #                                                               xmin = MLD.avg - MLD.sd,
    #                                                               xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
    #   }
    # }

    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP1, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }

    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP2, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }

    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP3, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }

    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP4, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }

    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP5, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }

    
    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP6, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }
    
    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP7, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }
 
    
    if("DEP" %in% params){
      avgPlot = avgPlot +
        geom_path(data = DEP8, aes(x = value, y =month, group = cluster),col="grey",linetype = 2)
      if(ribbon == TRUE){
        avgPlot = avgPlot + geom_ribbon(data = plotData[[2]], aes(x = MLD.avg, y = month,
                                                                  xmin = MLD.avg - MLD.sd,
                                                                  xmax = MLD.avg + MLD.sd, group=cluster), fill = "black", alpha = .08)
      }
    }

    ggsave(plot = avgPlot,
           paste("figures/",unique(plotData$PARAM)[i],outputname,".png",sep=""),
           width = 11.5, height = 3.5)
  }
}

