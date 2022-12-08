
ncImport_Arabian_B<-function(floatList = c("2902210.csv","2902211.csv","2902270.csv","2902271.csv","2902272.csv","2902275.csv","2902276.csv","2902277.csv"),
                   #Path = "Arabian/",
                   outputDirectory = "Arabian/raw_B/"){ 
  require(RNetCDF)
  require(zoo)
  require(curl)
  require(stringr)
  #require(tidyverse)
  require(lubridate)
  #ranger les flotteur, profils et indiquer ou ils sont (path)
  # floatData <- read_csv(floatList)%>% 
  #   group_by(directory,float,basin) %>% 
  #   summarize(n=n())
  
  for(i in 1:length(floatList)){
      print(paste("File",i,"of",length(floatList)))
      currentFloat <- floatList[i] %>% str_remove(".csv")
      currentBasin <- "Arabian"
      cycle <- list.files(paste("Arabian/B_files/",currentFloat,"/profiles",sep="")) %>% str_remove(paste("BR",currentFloat,"_",sep="")) %>% str_remove(".nc")
      
      for (j in 1:(length(cycle))){
          floatPath <- paste("Arabian/B_files/",currentFloat,"/profiles/BR",currentFloat,"_",cycle[j],".nc",sep="")
          
          argo <- read.nc(open.nc(floatPath)) #lire le fichier nc
          # check if BBP profile is good (>0)
          if(any(argo$BBP700 < 0,na.rm=T)){next}
          
          nRows<-nrow(argo$CHLA) #Rows correspond to number of depths for each profile
          nCols<-length(argo$CYCLE_NUMBER) #Cols correspond to total number of profiles
          
          #trouve les profils ou il y a tous ces parametres
          colsOfInterest=grep("^PRES|TEMP|PSAL|DOWN|DOXY|CHLA|BBP|CDOM|NITRATE",names(argo))
          colNames<-names(argo)[colsOfInterest]
          
          refDate<-argo$REFERENCE_DATE_TIME
          jDate<-argo$JULD
          jDateQC<-argo$JULD_QC
          jDate_location<-argo$JULD_LOCATION
          
          latitude<-rep(na.fill(argo$LATITUDE,"extend"),each = nRows)
          longitude<-rep(na.fill(argo$LONGITUDE,"extend"),each = nRows)
          
          cycleNumber = rep(argo$CYCLE_NUMBER,each=nRows)
          jDate<-rep(jDate,each=nRows)
          refDate<-rep(refDate,times=length(jDate))
          direction<-str_split(argo$DIRECTION,'')[[1]]
          
          argoSub<-argo[colsOfInterest]
          outputMatrix = matrix(nrow = nCols * nRows, ncol = length(names(argoSub)))
          
          for(k in 1:length(names(argoSub))){
            outputMatrix[,k] = vectorize(argoSub[[k]],colNames[k],nCols,nRows)
          }
          
          outputDF <- data.frame(outputMatrix,stringsAsFactors = FALSE)
          colnames(outputDF)=colNames
          
          outputDF<-cbind(basin=currentBasin,float=currentFloat,cycleNumber,longitude,latitude,
                          date = ddays(jDate)+ymd_hms(refDate),direction=rep(direction,each=nRows),
                          year = year(ddays(jDate)+ymd_hms(refDate)),outputDF)
          outputDF <- outputDF[(nRows+1):(nRows*2),]
          
          
          if(exists("totalDF")){
            totalDF <- rbind(totalDF,outputDF)
          }
          
          if(!exists("totalDF")){
            totalDF <- outputDF
          }
      }

      totalDF<-totalDF %>% 
        filter(direction=="A",
               as.numeric(PRES)<=1000) %>% 
        mutate(floatID = paste(float,year-min(year,na.rm=TRUE)+1,sep="_")) %>% 
        select(-direction)
      print("Writing CSV")
      write_csv(totalDF, paste(outputDirectory,currentFloat,".csv",sep=""))
      rm(totalDF)
  }
}

#################################################################################
#######################################################################################

ncImport_Arabian_S<-function(floatList = c("2902210.csv","2902211.csv","2902270.csv","2902271.csv","2902272.csv","2902275.csv","2902276.csv","2902277.csv"),
                           #Path = "Arabian/",
                           outputDirectory = "Arabian/raw_S/"){ 
  require(RNetCDF)
  require(zoo)
  require(curl)
  require(stringr)
  #require(tidyverse)
  require(lubridate)
  #ranger les flotteur, profils et indiquer ou ils sont (path)
  # floatData <- read_csv(floatList)%>% 
  #   group_by(directory,float,basin) %>% 
  #   summarize(n=n())
  
  for(i in 1:length(floatList)){
    print(paste("File",i,"of",length(floatList)))
    currentFloat <- floatList[i] %>% str_remove(".csv")
    currentBasin <- "Arabian"
    cycle <- list.files(paste("Arabian/S_files/",currentFloat,"/profiles",sep="")) %>% str_remove(paste("SR",currentFloat,"_",sep="")) %>% str_remove(".nc")
    
    for (j in 1:(length(cycle))){
      floatPath <- paste("Arabian/S_files/",currentFloat,"/profiles/SR",currentFloat,"_",cycle[j],".nc",sep="")
      
      argo <- read.nc(open.nc(floatPath)) #lire le fichier nc
      # check if BBP profile is good (>0)
      if(any(argo$BBP700 < 0,na.rm=T)){next}
      
      nRows<-length(argo$CHLA) #Rows correspond to number of depths for each profile
      nCols<-length(argo$CYCLE_NUMBER) #Cols correspond to total number of profiles
      
      #trouve les profils ou il y a tous ces parametres
      colsOfInterest=grep("^PRES|TEMP|PSAL|DOWN|DOXY|CHLA|BBP|CDOM|NITRATE",names(argo))
      colNames<-names(argo)[colsOfInterest]
      
      refDate<-argo$REFERENCE_DATE_TIME
      jDate<-argo$JULD
      jDateQC<-argo$JULD_QC
      jDate_location<-argo$JULD_LOCATION
      
      latitude<-rep(argo$LATITUDE,each = nRows)
      longitude<-rep(argo$LONGITUDE,each = nRows)
      
      cycleNumber = rep(argo$CYCLE_NUMBER,each=nRows)
      jDate<-rep(jDate,each=nRows)
      refDate<-rep(refDate,times=length(jDate))
      direction<-str_split(argo$DIRECTION,'')[[1]]
      
      argoSub<-argo[colsOfInterest]
      outputMatrix = matrix(nrow = nCols * nRows, ncol = length(names(argoSub)))
      
      for(k in 1:length(names(argoSub))){
        outputMatrix[,k] = vectorize(argoSub[[k]],colNames[k],nCols,nRows)
      }
      
      outputDF <- data.frame(outputMatrix,stringsAsFactors = FALSE)
      colnames(outputDF)=colNames
      
      outputDF<-cbind(basin=currentBasin,float=currentFloat,cycleNumber,longitude,latitude,
                      date = ddays(jDate)+ymd_hms(refDate),direction=rep(direction,each=nRows),
                      year = year(ddays(jDate)+ymd_hms(refDate)),outputDF)
      
      
      if(exists("totalDF")){
        totalDF <- rbind(totalDF,outputDF)
      }
      
      if(!exists("totalDF")){
        totalDF <- outputDF
      }
    }
    
    totalDF<-totalDF %>% 
      filter(direction=="A",
             as.numeric(PRES)<=1000) %>% 
      mutate(floatID = paste(float,year-min(year,na.rm=TRUE)+1,sep="_")) %>% 
      select(-direction)
    print("Writing CSV")
    write_csv(totalDF, paste(outputDirectory,currentFloat,".csv",sep=""))
    rm(totalDF)
  }
}



# file <- read_csv("data/rawCSV/2902086.csv")
# Arabian <- list.files("Arabian/raw") %>% str_remove(".csv")
# for (i in 1:8){
#   currentFloat <- Arabian[i]
#   file_Arabian <- read_csv(paste("Arabian/raw/",currentFloat,".csv",sep=""))  
#   file_Arabian <- file_Arabian[,which(colnames(file_Arabian) %in% colnames(file))]
#   write_csv(file_Arabian,paste("Arabian/raw/",currentFloat,".csv",sep=""))
# }
##################################################################################
floatList = c("2902210.csv","2902211.csv","2902270.csv","2902271.csv","2902272.csv","2902275.csv","2902276.csv","2902277.csv")

for(i in 1:length(floatList)){
  print(paste("File",i,"of",length(floatList)))
  currentFloat <- floatList[i] %>% str_remove(".csv")
  
  file_S <- read_csv(paste("Arabian/raw_S/",currentFloat,".csv",sep=""))
  file_B <- read_csv(paste("Arabian/raw_B/",currentFloat,".csv",sep=""))
  PRES_S <- file_S$PRES
  PRES_B <- file_B$PRES
  CHLA_S <- file_S$CHLA
  CHLA_B <- file_B$CHLA
  CHLA_ADJ_S <- file_S$CHLA_ADJUSTED
  CHLA_ADJ_B <- file_B$CHLA_ADJUSTED
  CHLA_QC_S <- file_S$CHLA_QC
  CHLA_ADJ_QC_S <- file_S$CHLA_ADJUSTED_QC
  CHLA_ADJ_QC_B <- file_B$CHLA_ADJUSTED_QC
  
  CHLA_ADJ_S[which(PRES_S %in% PRES_B & is.na(CHLA_S)==F)] <- CHLA_ADJ_B
  file_S$CHLA_ADJUSTED <- CHLA_ADJ_S
  
  i=1
  while(CHLA_B[i]==CHLA_S[which(PRES_S %in% PRES_B & is.na(CHLA_S)==F)][i]){
    i=i+1
  }
}


#################################################################################
csvQC_Arabian <- function (inputDirectory = "Arabian/raw/",
                           outputDirectory = "Arabian/qc/",
                           overwrite = FALSE,
                           unlink = FALSE){
  library(dplyr)
  library(lubridate)
  library(RNetCDF)
  library(stringr)
  library(readxl)
  
  files <- list.files(inputDirectory,pattern=".csv",full.names = T)
  
  for(i in 1:length(files)){
    print(paste("File",i,"of",length(files)))
    import<-read_csv(files[i],show_col_types = F) %>%
      mutate(date=ymd_hms(date))
    #View(import[import$cycleNumber==184,])
    
    currentFloat = import$float[1]
    cols = colnames(import)
    
    if(overwrite==FALSE){
      existingFiles <- list.files(outputDirectory)
      if(length(grep(currentFloat,existingFiles))!=0){
        next
      }
    }
    
    #   #ajuster les données de chla et bbp
    #   for(j in 1:length(cols)){
    #     print(j)
    #     paramColName = cols[j]
    #     #CHL_ADJ measurements are scaled by calibration parameters and should always be used over the CHLA values
    #     #this is somewhat different than other parameters, so there's a special case
    #     if(grepl("CHLA$", paramColName)){
    #       print(paramColName)
    #       pColNo = j
    #       qcColNo = grep(paste("^",paramColName,"_QC$",sep=""),cols)
    #       adjPColNo = grep(paste("^",paramColName,"_ADJUSTED$",sep=""),cols)
    #       adjQcColNo = grep(paste("^",paramColName,"_ADJUSTED_QC$",sep=""),cols)
    #       import[adjPColNo][import[adjQcColNo]==4]=NA #Set bad adjData to NA
    #       import[pColNo] = import[adjPColNo]
    #       import[qcColNo] = import[adjQcColNo]
    #       rm(pColNo,qcColNo,adjPColNo,adjQcColNo)
    #       next}
    #     #There are no BBP_ADJ measurements, so only BBP700 results should be used; another special case
    #     if(grepl("BBP700$", paramColName)){
    #       print(paramColName)
    #       pColNo = j
    #       qcColNo = grep(paste("^",paramColName,"_QC$",sep=""),cols)
    #       import[pColNo][import[qcColNo]==4 | is.na(import[qcColNo])==TRUE]=NA #Set bad data to NA
    #       rm(pColNo,qcColNo)
    #       next}
    #     if(paste(paramColName,"_QC",sep="") %in% cols &
    #        paste(paramColName,"_ADJUSTED",sep="") %in% cols &
    #        grepl("CHLA$", paramColName) == FALSE){
    #       print(paramColName)
    #       pColNo = j
    #       qcColNo = grep(paste("^",paramColName,"_QC$",sep=""),cols)
    #       adjPColNo = grep(paste("^",paramColName,"_ADJUSTED$",sep=""),cols)
    #       adjQcColNo = grep(paste("^",paramColName,"_ADJUSTED_QC$",sep=""),cols)
    #       import[pColNo][is.na(import[qcColNo]) == TRUE] = NA #Set param cols with PARAM_QC of NA to NA
    #       import[pColNo][import[qcColNo] == 3 | import[qcColNo] == 4] = NA #Set param cols with PARAM_QC of 3 or 4 to NA
    #       import[adjPColNo][import[qcColNo] == 1 | import[qcColNo] == 2] = NA #Set adjParam cols with PARAM_QC of 1 or 2 to NA
    #       import[adjPColNo][import[adjQcColNo]==4| is.na(import[adjQcColNo])==TRUE]=NA #Set adjParam cols with ADJ_PARAM_QC of 4 or NA to NA
    #       import[pColNo] = apply(import[,c(pColNo,adjPColNo)],1,mean,na.rm=TRUE)
    #       rm(pColNo,qcColNo,adjPColNo,adjQcColNo)
    #     }
    #     else{
    #       next}
    # }
    
    #ajuster les données de chla et bbp
    for(j in 1:length(cols)){
      print(j)
      paramColName = cols[j]
      #CHL_ADJ measurements are scaled by calibration parameters and should always be used over the CHLA values
      #this is somewhat different than other parameters, so there's a special case
      if(grepl("CHLA$", paramColName)){
        print(paramColName)
        pColNo = j
        qcColNo = grep(paste("^",paramColName,"_QC$",sep=""),cols)
        adjPColNo = grep(paste("^",paramColName,"_ADJUSTED$",sep=""),cols)
        adjQcColNo = grep(paste("^",paramColName,"_ADJUSTED_QC$",sep=""),cols)
        import[adjPColNo][import[adjQcColNo] == 3 | import[adjQcColNo] == 4]=NA #Set bad (3+4) adjData to NA
        import[pColNo] = import[adjPColNo] #deplacer la colonne adj dans la colonne 'normale'
        import[qcColNo] = import[adjQcColNo] #deplacer colonne qcadj dans colonne qc
        rm(pColNo,qcColNo,adjPColNo,adjQcColNo)
        next}
      #There are not BBP_ADJ measurements every time, so only BBP700 results should be used; another special case
      if(grepl("BBP700$", paramColName)){
        print(paramColName)
        pColNo = j
        qcColNo = grep(paste("^",paramColName,"_QC$",sep=""),cols)
        adjPColNo = grep(paste("^",paramColName,"_ADJUSTED$",sep=""),cols)
        adjQcColNo = grep(paste("^",paramColName,"_ADJUSTED_QC$",sep=""),cols)
        import[pColNo][import[qcColNo]==4 | is.na(import[qcColNo])==TRUE]=NA #Set bad data to NA pour le realtime
        import[adjPColNo][import[adjQcColNo] == 3 | import[adjQcColNo] == 4]=NA #Set bad (3+4) adjData to NA pour le adjusted
        if(is.na(import[adjPColNo])==FALSE){ #si il y a des paramètres BBP ajustés : les déplacer dans la colonne BBP
          import[pColNo] = import[adjPColNo]
          import[qcColNo] = import[adjQcColNo]
        }
        else{
          next}
        rm(pColNo,qcColNo,adjPColNo,adjQcColNo)
        next}
      
      if(paste(paramColName,"_QC",sep="") %in% cols &
         paste(paramColName,"_ADJUSTED",sep="") %in% cols &
         grepl("CHLA$", paramColName) == FALSE &
         grepl("BBP700$", paramColName) == FALSE){
        print(paramColName)
        pColNo = j
        qcColNo = grep(paste("^",paramColName,"_QC$",sep=""),cols)
        adjPColNo = grep(paste("^",paramColName,"_ADJUSTED$",sep=""),cols)
        adjQcColNo = grep(paste("^",paramColName,"_ADJUSTED_QC$",sep=""),cols)
        import[pColNo][is.na(import[qcColNo]) == TRUE] = NA #Set param cols with PARAM_QC of NA to NA
        import[pColNo][import[qcColNo] == 3 | import[qcColNo] == 4] = NA #Set param cols with PARAM_QC of 3 or 4 to NA
        import[adjPColNo][import[qcColNo] == 1 | import[qcColNo] == 2] = NA #Set adjParam cols with PARAM_QC of 1 or 2 to NA
        import[adjPColNo][import[adjQcColNo]==4| is.na(import[adjQcColNo])==TRUE]=NA #Set adjParam cols with ADJ_PARAM_QC of 4 or NA to NA
        import[pColNo] = apply(import[,c(pColNo,adjPColNo)],1,mean,na.rm=TRUE)
        rm(pColNo,qcColNo,adjPColNo,adjQcColNo)
      }
      else{
        next}
    }
    
    import=import[-grep("_QC$|_ADJUSTED$|_ERROR$",cols)] #Remove QC columns
    
    ### %m+% ensures that number of days in a month aren't exceeded when adding 6 months (as in 
    ### 2004-01-31 + months(1))
    
    ###Create column of shifted dates; dates for Southern Hemisphere profiles will be shifted 6 months;
    ###dates for Northern Hemisphere profiles are not be changed
    ###This only changes dates for floats that only occupy the Southern Hemisphere. Floats crossing the 
    ###equator will not be affected
    import$shift = NA
    import$jDay = yday(import$date)
    import$shift = 0
  
    
    condensed <- import %>% 
      group_by(cycleNumber) %>% 
      #There are 1 or 2 floats with completely missing PSAL/TEMP data; the filter
      #settings remove all corresponding data
      filter(all(is.na(date))==FALSE,
             all(is.na(PRES))==FALSE,
             all(is.na(PSAL))==FALSE,
             all(is.na(TEMP))==FALSE,
             all(is.na(CHLA))==FALSE,
             all(is.na(BBP700))==FALSE)
    
    #enlever les profils avec moins de 5 lignes
    if(nrow(condensed)<5){
      print("Reject!")
      next}
    
    final <- condensed %>% 
      filter(is.na(PRES)==FALSE) %>%
      group_by(date) %>% 
      mutate(minPres = min(PRES,na.rm=TRUE),
             maxPres = max(PRES,na.rm=TRUE),
             depths = length(which(PRES < 300)))%>% 
      filter(depths>15,
             minPres < 5,
             maxPres > 300)
    
    
    if(nrow(final)<5){
      print("Reject!")
      next}
    
    #Function from Stack exchange https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame/18143097#18143097
    #Replaces NaN cells introudced while averaging multiple NA values in data adjustment section to NA values
    is.nan.data.frame <- function(x){
      do.call(cbind, lapply(x, is.nan))}
    
    final[is.nan.data.frame(final)] <- NA
    
    write_csv(final,paste(outputDirectory,currentFloat,".csv",sep=""))
    rm(import,condensed,final)
    if(unlink==TRUE){
      unlink(files[i])
    }
  }
}

###############################################################################

csvAdjust_Arabian <- function(inputDirectory = "data/qcCSV/",
                      outputDirectory = "data/finalCSV/",
                      overwrite = FALSE,
                      unlink = FALSE) {
  
  #require(tidyverse)
  require(tidyr)
  require(lubridate)
  require(readxl)
  require(cowplot)
  require(gsw)
  require(rworldmap)
  require(zoo)
  require(oce)
  require(GCalignR)
  
  files <- list.files(inputDirectory,pattern=".csv",full.names = T)
  
  for(i in 1:length(files)){
    print(paste("File",i,"of",length(files)))
    import = read_csv(files[i])
    currentBasin = import$basin[1]
    currentFloat = import$float[1]
    
    if(overwrite==FALSE){
      existingFiles <- list.files(outputDirectory)
      if(length(grep(currentFloat,existingFiles))!=0){
        next
      }
    }
    
    #Identifies parameter columns
    paramCols <- which(colnames(import) %in% c("CHLA","BBP700","BBP532","DOXY","NITRATE","PSAL","TEMP",
                                               "DOWNWELLING_PAR","DOWN_IRRADIANCE380",
                                               "DOWN_IRRADIANCE412",
                                               "DOWN_IRRADIANCE490","CDOM"))
    
    
    #Identifies depth adjustment columns                                          
    #dPresCols <- which(colnames(import) %in% paste(colnames(import)[paramCols],"_dPRES",sep=""))
    
    #Creates a dataframe containing vectors for dates, depths, adjustement factors,
    #parameter names, and parameter values, 
    dateVector <- rep(import$date,times = length(paramCols))
    yearVector <- rep(year(import$date),times = length(paramCols))
    presVector <- rep(import$PRES,times = length(paramCols))
    paramNameVector <- rep(colnames(import)[paramCols],each = nrow(import))
    paramVector <- as.vector(unlist(import[,c(paramCols)]))
    #dPresVector <- as.vector(unlist(import[,c(dPresCols)]))
    
    metaData <- import[c(-paramCols)] %>% 
      group_by(floatID,date,basin) %>% 
      summarize_all(mean) %>% 
      select(-PRES)
    
    #Merges above vectors into a single data frame; removes missing parameter values
    #and bins by depth
    originalData <- data.frame(date = dateVector, year = yearVector, PRES = presVector,
                               param = paramNameVector,
                               value = paramVector,row.names = NULL) %>% 
      mutate(timeSeries = year - min(year)+1,
             PRES = PRES) %>% 
      group_by(date,param) %>% 
      #Removes all NA values; NOTE! If there are no measurements for any param at 
      #less than 5 meters, this can create profiles with minimum depths of more than 5m
      #this can conflict with parameter estimates in the csvFinal script
      filter(is.na(value)==FALSE) %>% 
      mutate(PRES = ceiling(PRES),
             n = sum(PRES<=300),
             unique = length(unique(value[PRES<=300])))%>% 
      #Removes any profile with fewer than 15 measurements remaining after removing NA
      filter(n > 15) %>% 
      group_by(date) %>% 
      #Recalculating minPres removes profiles lacking measurements at less than 5m;
      #this is redundant, given that these profiles should have been removed in 
      #bulkCsvClean, but...
      mutate(minPres = min(PRES,na.rm=TRUE)) %>% 
      filter(minPres <= 10) %>% 
      select(-minPres) %>% 
      group_by(timeSeries,date,PRES,param) %>% 
      #Averages by depth
      summarize(value = mean(value,na.rm=TRUE)) %>% 
      ungroup()
    
    adjData<-originalData
    originalData$param <- paste(originalData$param,"_ORIG",sep="")
    originalData <- originalData %>% 
      spread(param,value) %>% 
      select(-timeSeries)
    
    #Creates subset of BBP values if BBP is amongst measured parameters
    #This is necessary to apply different functions to BBP and other columns
    if("BBP700" %in% colnames(import)){
      nonBBP <- adjData %>% 
        filter(param != "BBP700") %>% 
        group_by(date,param) %>% 
        mutate(value = rollmedian(value,3,fill = "extend"),
               value = rollmean(value,3,fill="extend"))
      
      bbpVector <- adjData %>% 
        filter(param == "BBP700") %>% 
        group_by(timeSeries) %>% 
        mutate(bbpThresh = quantile(value,0.995)) %>% 
        filter(value < bbpThresh) %>% 
        group_by(date,param) %>% 
        mutate(value = rollmedian(value,7,fill = "extend"),
               value = rollmean(value,5,fill = "extend")) %>% 
        select(-bbpThresh)
      
      adjData <- rbind(bbpVector,nonBBP) %>% 
        spread(param,value) %>% 
        ungroup()
    }
    #If BBP is not amongst measured parameters, just spreads the existing dataVector  
    if("BBP700" %in% colnames(import)==FALSE){
      adjData<-adjData %>% 
        group_by(date,param) %>% 
        mutate(value = rollmedian(value,3,fill = "extend"),
               value = rollmean(value,3,fill="extend")) %>% 
        spread(param,value) %>% 
        ungroup()
    }
    
    #This loop interpolates each profile for each parameter; two loops are required
    #in order to catch profiles where all values are missing (note that profiles
    #with fewer than 15 measurements are set to NA when creating dataVector)
    for(i in unique(adjData$date)){
      dateRange = which(adjData$date == i)
      for(j in 4:ncol(adjData)){
        if(all(is.na(adjData[dateRange,j]))){
          next}
        else{
          adjData[dateRange,j] = na.fill(adjData[dateRange,j], fill = "extend")
        }
      }
    }
    
    
    adjData<-as.data.frame(adjData)
    originalData<-as.data.frame(originalData)
    adjOrig <- merge(originalData,adjData,by=c("date","PRES"))
    
    #Merges interpolated data with import file
    adjImport<-metaData %>% 
      merge(adjOrig,by=c("date"),all.y=TRUE) %>% 
      group_by(date) %>% 
      #Following lines fill in any empty cells created by merge
      mutate(float = na.fill(float,fill="extend"),
             year = na.fill(year,fill="extend"),
             #jDay = na.fill(jDay,fill="extend"),
             latitude = na.fill(latitude,fill="extend"),
             longitude = na.fill(longitude,fill="extend"),
             SIG_ORIG = gsw_sigma0(PSAL_ORIG,TEMP_ORIG),
             SIG = gsw_sigma0(PSAL,TEMP)) %>% 
      group_by(float) %>% 
      mutate(floatID = paste(float,"_", timeSeries,sep=""))
    
    adjImport<-data.frame(adjImport)
    rm(adjOrig,originalData,adjData,import)
    adjImport$jDay = yday(adjImport$date)
    #adjImport$density = swRho(adjImport$PSAL,adjImport$TEMP,0) #potential density
    # adjImport$MLD = NA #MLD depth
    # for(j in 1:max(adjImport$cycleNumber)){
    #   if(j %in% adjImport$cycleNumber){
    #     adjImport$MLD[which(adjImport$cycleNumber==j)] <- mld(x = adjImport$density[which(adjImport$cycleNumber==j)],
    #                                                     depth = adjImport$PRES[which(adjImport$cycleNumber==j)],
    #                                                     ref.depths = 10,
    #                                                     criteria = 0.03,
    #                                                     default.depth=NA)
    #   }
    # }
    #if("DOXY" %in% colnames(adjImport)){
    #adjImport$AOU = oxySat(adjImport$DOXY,adjImport$TEMP,adjImport$PSAL)} #AOU
    write_csv(adjImport,
              paste(outputDirectory,currentFloat,".csv",sep=""))
    if(unlink==TRUE){
      unlink(files[i])
    }
  }
}

##############################################################################
csvAdjust <- function(inputDirectory = "data/qcCSV/",
                      outputDirectory = "data/finalCSV/",
                      overwrite = FALSE,
                      unlink = FALSE) {
  
  #require(tidyverse)
  require(tidyr)
  require(lubridate)
  require(readxl)
  require(cowplot)
  require(gsw)
  require(rworldmap)
  require(zoo)
  require(oce)
  require(GCalignR)
  
  files <- list.files(inputDirectory,pattern=".csv",full.names = T)
  
  for(i in 1:length(files)){
    print(paste("File",i,"of",length(files)))
    import = read_csv(files[i],show_col_types = F)
    currentBasin = import$basin[1]
    currentFloat = import$float[1]
    
    if(overwrite==FALSE){
      existingFiles <- list.files(outputDirectory)
      if(length(grep(currentFloat,existingFiles))!=0){
        next
      }
    }
    
    #Identifies parameter columns
    paramCols <- which(colnames(import) %in% c("CHLA","BBP700","BBP532","DOXY","NITRATE","PSAL","TEMP",
                                               "DOWNWELLING_PAR","DOWN_IRRADIANCE380",
                                               "DOWN_IRRADIANCE412",
                                               "DOWN_IRRADIANCE490","CDOM"))
    
    
    #Identifies depth adjustment columns                                          
    dPresCols <- which(colnames(import) %in% paste(colnames(import)[paramCols],"_dPRES",sep=""))
    
    #Creates a dataframe containing vectors for dates, depths, adjustement factors,
    #parameter names, and parameter values, 
    dateVector <- rep(import$date,times = length(paramCols))
    yearVector <- rep(year(import$date),times = length(paramCols))
    presVector <- rep(import$PRES,times = length(paramCols))
    paramNameVector <- rep(colnames(import)[paramCols],each = nrow(import))
    paramVector <- as.vector(unlist(import[,c(paramCols)]))
    dPresVector <- as.vector(unlist(import[,c(dPresCols)]))
    
    metaData <- import[c(-paramCols,-dPresCols)] %>% 
      group_by(floatID,date,basin) %>% 
      summarize_all(mean) %>% 
      select(-PRES)
    
    #Merges above vectors into a single data frame; removes missing parameter values
    #and bins by depth
    originalData <- data.frame(date = dateVector, year = yearVector, PRES = presVector,
                               param = paramNameVector,
                               value = paramVector, dPres = dPresVector,row.names = NULL) %>% 
      mutate(timeSeries = year - min(year)+1,
             PRES = PRES+dPres) %>% 
      group_by(date,param) %>% 
      #Removes all NA values; NOTE! If there are no measurements for any param at 
      #less than 5 meters, this can create profiles with minimum depths of more than 5m
      #this can conflict with parameter estimates in the csvFinal script
      filter(is.na(value)==FALSE) %>% 
      mutate(PRES = ceiling(PRES),
             n = sum(PRES<=300),
             unique = length(unique(value[PRES<=300])))%>% 
      #Removes any profile with fewer than 15 measurements remaining after removing NA
      filter(n > 15) %>% 
      group_by(date) %>% 
      #Recalculating minPres removes profiles lacking measurements at less than 5m;
      #this is redundant, given that these profiles should have been removed in 
      #bulkCsvClean, but...
      mutate(minPres = min(PRES,na.rm=TRUE)) %>% 
      filter(minPres <= 5) %>% 
      select(-minPres) %>% 
      group_by(timeSeries,date,PRES,param) %>% 
      #Averages by depth
      summarize(value = mean(value,na.rm=TRUE)) %>% 
      ungroup()
    
    adjData<-originalData
    originalData$param <- paste(originalData$param,"_ORIG",sep="")
    originalData <- originalData %>% 
      spread(param,value) %>% 
      select(-timeSeries)
    
    #Creates subset of BBP values if BBP is amongst measured parameters
    #This is necessary to apply different functions to BBP and other columns
    if("BBP700" %in% colnames(import)){
      nonBBP <- adjData %>% 
        filter(param != "BBP700") %>% 
        group_by(date,param) %>% 
        mutate(value = rollmedian(value,3,fill = "extend"),
               value = rollmean(value,3,fill="extend"))
      
      bbpVector <- adjData %>% 
        filter(param == "BBP700") %>% 
        group_by(timeSeries) %>% 
        mutate(bbpThresh = quantile(value,0.995)) %>% 
        filter(value < bbpThresh) %>% 
        group_by(date,param) %>% 
        mutate(value = rollmedian(value,7,fill = "extend"),
               value = rollmean(value,5,fill = "extend")) %>% 
        select(-bbpThresh)
      
      adjData <- rbind(bbpVector,nonBBP) %>% 
        spread(param,value) %>% 
        ungroup()
    }
    #If BBP is not amongst measured parameters, just spreads the existing dataVector  
    if("BBP700" %in% colnames(import)==FALSE){
      adjData<-adjData %>% 
        group_by(date,param) %>% 
        mutate(value = rollmedian(value,3,fill = "extend"),
               value = rollmean(value,3,fill="extend")) %>% 
        spread(param,value) %>% 
        ungroup()
    }
    
    #This loop interpolates each profile for each parameter; two loops are required
    #in order to catch profiles where all values are missing (note that profiles
    #with fewer than 15 measurements are set to NA when creating dataVector)
    for(i in unique(adjData$date)){
      dateRange = which(adjData$date == i)
      for(j in 4:ncol(adjData)){
        if(all(is.na(adjData[dateRange,j]))){
          next}
        else{
          adjData[dateRange,j] = na.fill(adjData[dateRange,j], fill = "extend")
        }
      }
    }
    
    
    adjData<-as.data.frame(adjData)
    originalData<-as.data.frame(originalData)
    adjOrig <- merge(originalData,adjData,by=c("date","PRES"))
    
    #Merges interpolated data with import file
    adjImport<-metaData %>% 
      merge(adjOrig,by=c("date"),all.y=TRUE) %>% 
      group_by(date) %>% 
      #Following lines fill in any empty cells created by merge
      mutate(float = na.fill(float,fill="extend"),
             year = na.fill(year,fill="extend"),
             jDay = na.fill(jDay,fill="extend"),
             latitude = na.fill(latitude,fill="extend"),
             longitude = na.fill(longitude,fill="extend"),
             SIG_ORIG = gsw_sigma0(PSAL_ORIG,TEMP_ORIG),
             SIG = gsw_sigma0(PSAL,TEMP)) %>% 
      group_by(float) %>% 
      mutate(floatID = paste(float,"_", timeSeries,sep=""))
    
    adjImport<-data.frame(adjImport)
    rm(adjOrig,originalData,adjData,import)
    adjImport$jDay = yday(adjImport$date)
    #adjImport$POC <- 9.776*10^4*(adjImport$BBP700)^1.166
    #if("DOXY" %in% colnames(adjImport)){
    #adjImport$AOU = oxySat(adjImport$DOXY,adjImport$TEMP,adjImport$PSAL)} #AOU
    write_csv(adjImport,
              paste(outputDirectory,currentFloat,".csv",sep=""))
    if(unlink==TRUE){
      unlink(files[i])
    }
  }
}

################################################################################
mld<-function(fileList,plot=FALSE,
              inputDirectory = "data/finalCSV/",
              outputDirectory = "data/finalCSV/"){
  require(gsw)
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read_csv(paste(inputDirectory,fileList[i],sep=""),show_col_types = F)
    
    temp$MLD_003 = NA
    for(j in 1:length(unique(temp$date))){
      print(j)
      
      range = which(temp$date==unique(temp$date)[j])
      if(length(range)==1){
        next
      }
      currentFloat = temp$floatID[range][1]
      temp$MLD_003[range] = sigma_mixedDepth(sal = temp$PSAL[range],
                                             temp = temp$TEMP[range],
                                             dbar = temp$PRES[range],
                                             thresh = 0.03,
                                             ref = 5)
      
      if(plot==TRUE){
        output = ggplot() +
          geom_line(aes(x = temp$PRES[range],y = temp$SIG[range])) +
          geom_vline(aes(xintercept = temp$MLD_003[range][1]), col = "red") +
          coord_flip() +
          scale_x_reverse() +
          theme_bw()
        
        ggsave(plot=output,
               paste("figures/parameterTestPlots/mldTest/",
                     currentFloat,"_",j,".png",sep=""))
      }
    }
    write_csv(temp, paste(outputDirectory,fileList[i],sep=""))
  }
}
##############################################################################
sigma_mixedDepth <- function(sal,temp,dbar,thresh,ref) {
  require(gsw)
  #browser()
  #Do not calculate if dbar starts at a depth greater than 12m, if measurements are taken at intervals >30m,
  #or if fewer than 5 measurements are available
  #Reference depth is typically 10 m
  
  if(all(is.na(sal))|all(is.na(temp))|all(is.na(dbar))){
    return(NA)
  }
  
  depthTemp <- approx(dbar,sal,1:max(dbar,na.rm=TRUE))$x
  #Use NA fill to replace missing values in cases where first measurement was taken
  #from a depth greater than 1 meter
  salTemp <-  na.fill(approx(dbar,sal,1:max(dbar,na.rm=TRUE))$y,"extend")
  tempTemp <-  na.fill(approx(dbar,temp,1:max(dbar,na.rm=TRUE))$y,"extend")
  
  sigma <- gsw_sigma0(salTemp,tempTemp)
  
  #Values are interpolated to 1m, so sigmaRef is just sigma at the reference value
  sigmaRef <- sigma[ref] #Density at reference depth
  #MLD is calculated as the depth at which density deviates from the reference density by 0.03 kg m^-3. 
  deltaSigma <- sigma - sigmaRef #This will equal 0 if there is a measurement exactly at the threshold
  #all negative values will be below the threshold; so the MLD is the smallest positive result of this equation
  MLD = depthTemp[deltaSigma>=0.03][1]
  return(MLD)
}
###############################################################################
province<-function(fileList,
                   inputDirectory = "data/finalCSV/",
                   outputDirectory = "data/finalCSV/"){
  for(i in 1:length(fileList)){
    print(paste(i,fileList[i]))
    temp<-read_csv(paste(inputDirectory,fileList[i],sep=""),show_col_types = F)
    
    currentFloat<-temp$float[1]
    
    temp<-temp%>% 
      group_by(date) %>% 
      mutate(province = zone(lat = mean(latitude,na.rm=TRUE),
                             lon = mean(longitude,na.rm=TRUE),
                             temp = TEMP_ORIG,
                             dep_temp = PRES))
    
    #if any province is NA : define the province with the one of the previous cycle NUmber
    if(any(is.na(temp$province)) & !all(is.na(temp$province)) & !is.na(temp$province[1])){
      NAs <- unique(temp$cycleNumber[which(is.na(temp$province))])
      cycles <- unique(temp$cycleNumber)
      ref_cycles <- which(cycles %in% NAs)-1
      for(j in 1:length(NAs)){
        temp$province[which(temp$cycleNumber==NAs[j])] <- unique(temp$province[(which(temp$cycleNumber==cycles[ref_cycles[j]]))])
      }
    }
    
    if(all(is.na(temp$province)) | is.na(temp$province[1])){
      temp$province[which(is.na(temp$province))] <- "to be determined"
      print("zone to be determined")
    }
    
    write_csv(temp,
              paste(outputDirectory,fileList[i],sep=""))
  }
}
###############################################################################
zone<-function(lat,lon,temp,dep_temp) {
  #browser()
  zone<-NA
  
  temp <- tibble(prof = dep_temp, temp = temp) %>% arrange(prof) %>% pull(temp)
  dep_temp <- tibble(prof = dep_temp, temp = temp) %>% arrange(prof) %>% pull(prof)
  
  if(all(is.na(temp))){
    zone <- NA
    return(zone)
  }
  
  if(lat > 45 & lon < -100  ) {
    zone<-"NPSPG"
  }
  
  # 23/02/19 : replaced below
  # if(lat > 50 & lat < 65 & lon > -65  & lon < -12.5  ) {
  #   zone<-"NASPG"
  # }
  
  if(lat > 60 & lat < 80 & lon > -12.5  & lon < 20  ) {
    zone<-"NS"
  }
  
  if(lat > 17.5 & lat < 32.5 & lon > -75  & lon < -20  ) {
    zone<-"NASTG"
  }
  
  if(lat > 17.5 & lat < 32.5 & lon > -175  & lon < -140  ) {
    zone<-"NPSTG"
  }
  
  
  #26/02/19 change lat max from -15 to -14 & lon max from -120 to -90
  if(lat > -30 & lat < -14 & lon > -180  & lon < -90  ) {
    zone<-"SPSTG"
  }
  
  #07/03/19 exlude Marquises (change lon min from -150 to -125)
  #26/02/19 add UPW
  if(lat > -14 & lat < 6 & lon > -125  & lon < -75  ) {
    zone<-"UPW"
  }
  
  if(lat > 0 & lat < 17.5 & lon > -35  & lon < -10  ) {
    zone<-"ASEW"
  }
  
  if(lat > 0 & lat < 17.5 & lon > -140  & lon < -100  ) {
    zone<-"PSEW"
  }
  
  if(lat > -30 & lat < -15 & lon > -45  & lon < 0  ) {
    zone<-"SASTG"
  }
  
  #26/02/19 Change lat max from -10 to -15
  #26/02/19 Change lon max from 100 to 110
  if(lat > -30 & lat < -15 & lon > 50  & lon < 110  ) {
    zone<-"SISTG"
  }
  
  if(lat > -5 & lat < 7.5 & lon > 40  & lon < 100  ) {
    zone<-"IEQ"
  }
  
  #27/02/19 Change lon min from 65 to 50  
  #26/02/19 ADD MOONS
  if(lat > -15 & lat < -5 & lon > 50  & lon < 110  ) {
    zone<-"MOONS"
  }
  
  # 26/02/19 Change lat min from > to >=
  if(lat >= 7.5 & lat < 24 & lon > 43  & lon < 100  ) {
    zone<-"IOMZ"
  }
  
  if(lat > 35 & lat < 43 & lon > -5  & lon < 6  ) {
    zone<-"WMS"
  }
  
  # 26/02/19 change lon min from 5 to 6
  if(lat > 35 & lat < 45 & lon > 5  & lon < 16 ) {
    zone<-"WMS"
  }
  
  if(lat > 30 & lat < 45 & lon > 16  & lon < 22.5  ) {
    zone<-"EMS"
  }
  
  if(lat > 30 & lat < 40 & lon > 22.5  & lon < 35  ) {
    zone<-"EMS"
  }
  
  if(lat > 41 & lat < 46 & lon > 12  & lon < 15.5  ) {
    zone<-"EMS"
  }
  
  if(lat > 30 & lat < 37.5 & lon > 10  & lon < 15.5  ) {
    zone<-"EMS"
  }
  
  #23/02/19
  if(lat > 40.98 & lat < 47.35 & lon > 27.67  & lon < 42.5  ) {
    zone<-"BKS"
  }
  
  #23/02/19
  if(lat > 10.78 & lat < 30.86 & lon > 32.78  & lon < 42.4  ) {
    zone<-"RDS"
  }
  
  #26/02/19 Add AUS
  if(lat > -30 & lat < -7 & lon > 110  & lon < 160  ) {
    zone<-"AUS"
  }
  
  #14/11/19 Modofy lon max from 125 to 120
  #26/02/19 Add CHINSS
  if(lat > 0 & lat < 24 & lon > 100  & lon < 120  ) {
    zone<-"CHINS"
  }
  
  #26/02/19 Add ARCH
  if(lat > -30 & lat < -5 & lon > 142  & lon < 180  ) {
    zone<-"ARCH"
  }
  
  #27/02/19 Convert lon min from -130 to 135
  #26/02/19 Add CALCUR
  if(lat > 30 & lat < 50 & lon > -135  & lon < -115  ) {
    zone<-"CAL"
  }
  
  if(lat > 0 & lat < 17.5 & lon > -140  & lon < -100  ) {
    zone<-"PSEW"
  }
  
  #23/02/19 
  if (lat < 65 & lat > 30 & lon < -10 & lon > -80) { # North Atlantic delimitation
    
    if(length(dep_temp) > 5 & 
       max(dep_temp,na.rm=T) >= 100 & 
       min(dep_temp,na.rm=T) <= 100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      #NASTG
      if (t_100 > 17.5) {
        zone<-"NASTG"
      }
      
      #NAC
      if (t_100 < 17.5 & t_100 > 8 ) {
        zone<-"NAC"
      }
      
      #NASPG
      if (t_100 < 8 ) {
        zone<-"NASPG"
      }
    }
  }
  
  #26/02/19
  if (lat > 65 ) { # Arctic delimitation
    
    if(length(dep_temp) > 5 & 
       max(dep_temp,na.rm=T) >= 100 & 
       min(dep_temp,na.rm=T) <= 100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      #ARCT
      if (t_100 < 3) {
        zone<-"ARCT"
      }
    }
  }
  
  
  # 14/11/19 Modify lon min from 130 to 120
  # 14/11/19 Modify lat min from 25 to 0
  #26/02/19 
  if (lat < 65 & lat > 0 & lon > 120) { # 
    
    if(length(dep_temp) > 5 & 
       max(dep_temp,na.rm=T) >= 100 & 
       min(dep_temp,na.rm=T) <= 100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      #NPSTG
      if (t_100 > 15) {
        zone<-"NPSTG"
      }
      
      #KURIOSHO
      if (t_100 < 15 & t_100 > 7.5 ) {
        zone<-"KURIO"
      }
      
      #NPSPG
      if (t_100 < 7.5 ) {
        zone<-"NPSPG"
      }
    }
  }
  
  
  if(lat > 65 & lat < 70 & lon > -72  & lon < -50  ) {
    zone<-"BAFF"
  }
  
  if(lat > 70 & lat < 80 & lon > -80  & lon < -50  ) {
    zone<-"BAFF"
  }
  
  # 26/02/19 Convert close to South Africa to STZ
  if(is.na(zone) &  lat > -30 & lat < -25 & lon > 0  & lon < 50  ) {
    zone<-"STZ"
  }
  
  
  if(lat < - 30) {
    
    if (dep_temp[max(which(!is.na(temp)))] >= 400 & min(dep_temp,na.rm=T)<=100) {
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      t_400<-NA
      t_400<-approx(dep_temp,temp,400)$y
      
      t_0_200<-NA
      t_0_200<-min(temp[which(dep_temp<=200)],na.rm=T)
      
      #STZ
      if (t_100 >= 11) {
        zone<-"STZ"
      }
      
      #SAZ
      if (t_100 < 11 & t_400 >= 5) {
        zone<-"SAZ"
      }
      
      #PFZ
      if (t_400 < 5 & t_0_200 >= 2) {
        zone<-"PFZ"
      }
      
      #ASZ_SIZ
      if (t_0_200 < 2) {
        zone<-"ASZ_SIZ"
      }
    }
    
    if (dep_temp[max(which(!is.na(temp)))] >= 380 & dep_temp[max(which(!is.na(temp)))] < 400 & min(dep_temp,na.rm=T)<=100) { #si derniere valeur de T entre 380 et 400m alors la prednre comme valeur pour 400
      t_100<-NA
      t_100<-approx(dep_temp,temp,100)$y
      
      t_400<-NA
      t_400<-temp[max(which(!is.na(temp)))]
      
      t_0_200<-NA
      t_0_200<-min(temp[which(dep_temp<=200)],na.rm=T)
      
      #STZ
      if (t_100 >= 11) {
        zone<-"STZ"
      }
      
      #SAZ
      if (t_100 < 11 & t_400 >= 5) {
        zone<-"SAZ"
      }
      
      #PFZ
      if (t_400 < 5 & t_0_200 >= 2) {
        zone<-"PFZ"
      }
      
      #ASZ_SIZ
      if (t_0_200 < 2) {
        zone<-"ASZ_SIZ"
      }
    }
  }
  
  #26/02/19 change lat max from -15 to -14 & lon max from -120 to -90
  if(lat > -30 & lat < -14 & lon > -180  & lon < -90  ) {
    zone<-"SPSTG"
  }
  
  return(zone)
}

############################################################################


###############################################################################"
#####################################LAUNCH####################################
##################################################################################
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

ncImport_Arabian(floatList = c("2902210.csv","2902211.csv","2902270.csv","2902271.csv","2902272.csv","2902275.csv","2902276.csv","2902277.csv"),
         outputDirectory = "Arabian/",
         overwrite=TRUE)

csvQC_Arabian(inputDirectory = "Arabian/raw_S/",
              outputDirectory = "Arabian/qc/",
              overwrite = FALSE,
              unlink = FALSE)

csvAdjust(inputDirectory = "Arabian/qc/",
          outputDirectory = "Arabian/final/",
          overwrite=TRUE,unlink=FALSE)

mld(list.files("Arabian/final",pattern=".csv",full.names = F),plot=F,inputDirectory = "Arabian/final/",outputDirectory = "Arabian/final/")
province(fileList=list.files("Arabian/final",pattern=".csv",full.names = F),inputDirectory = "Arabian/final/",outputDirectory = "Arabian/") #charge zone and par_prof functions
parEstimator(inputDirectory = "Arabian/", outputDirectory = "Arabian/parData/", plot=F)
isolumes(fileList<-list.files("Arabian/final/",full.names = F,pattern=".csv"),
         parFiles<-list.files("Arabian/parData",pattern=".csv",full.names = T),
         inputDirectory = "Arabian/", 
         outputDirectory = "Arabian/")
