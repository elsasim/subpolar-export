################################################################################
################################################################################
################################################################################
getMetaData<-function(params = c("CHLA","BBP700"),
                      includedBasins = NA,
                      excludedBasins = NA, 
                      maxElevation = NA,
                      timeseries = TRUE){
  require(curl)
  require(lubridate)
  require(stringr)
  #require(tidyverse)
  require(dplyr)
  
  #Create temporary file as tmp
  floatData <- fread("/home/master/Documents/index/argo_bio-profile_index.txt")
  
  file_split <- str_split(floatData$file,"/",simplify=T)
  floatData$directory <- file_split[,1]
  floatData$float <- file_split[,2]
  floatData$date_update <- ymd_hms(floatData$date_update)
  floatData$date<-ymd_hms(floatData$date)
  
  #Check columns for presence/absence of PRES, CHLA, BBP; add column for duration between
  #first and last profile for a given float
  for(i in 1:length(params)){
    temp<-data.frame(grepl(params[i],floatData$parameters))
    colnames(temp) = params[i]
    floatData<-cbind(floatData,temp)
    floatData<-floatData[-which(temp==FALSE),]
    file_split<-file_split[-which(temp==FALSE),]
  }
  
  floatData$year = year(floatData$date)
  floatData$jDay = yday(floatData$date)
  floatData<-floatData %>% 
    mutate(timeseries = year - min(year,na.rm = TRUE),
           floatID = paste(float,timeseries,sep=""))
  
  params_split <- str_split(floatData$parameters," ",simplify=T) #spliter parametres
  data_mode_split <- str_split(floatData$parameter_data_mode,"",simplify=T) #spliter data_mode
  
  NoCol <- vector()
  bad_CHLA <- vector()
  
  #selectionner les profils de CHLA qui ne sont pas en adjusted ou delayed :
  for(i in 1:nrow(floatData)){
    param_CHLA = "CHLA" %in% params_split[i,]
    if(param_CHLA) {
      NoCol[i] <- which(params_split[i,] == "CHLA")
      if (data_mode_split[i,NoCol[i]]=="D"|data_mode_split[i,NoCol[i]]=="A"){
        bad_CHLA[i] = NA
      }
      if(data_mode_split[i,NoCol[i]]!="D" & data_mode_split[i,NoCol[i]]!="A"){
        bad_CHLA[i] <- i
      }
    }
    if(param_CHLA == FALSE) {
    }
  }
  bad_CHLA <- bad_CHLA[!is.na(bad_CHLA)]
  floatData <- floatData[-bad_CHLA,]
  # Il reste les profils avec CHLA = A/D et quand pas CHLA
  
  floatData <- floatData[-which(is.na(floatData$date)),]
  
  if(timeseries==TRUE){
    print("Defining basins...")
    floatData<-floatData %>% 
      group_by(float) %>% 
      mutate(basin = basinFinder(median(longitude,na.rm=TRUE),median(latitude,na.rm=TRUE))) %>% 
      group_by(floatID) %>% 
      mutate(profiles = n(),
             minJday = min(jDay,na.rm=TRUE),
             maxJday = max(jDay,na.rm=TRUE)) %>% 
      #Only keep floats with more than 32 profiles (the minimum required to provide 
      #an annual time series with sampling resolution of 10 days, allowing for 2 missing
      #profiles at the start or end of the year)
      filter(profiles >= 32,
             minJday <= 25,
             maxJday >= 340)
    
    if(is.na(includedBasins)==FALSE){
      floatData <- floatData[which(floatData$basin %in% includedBasins),]
    }
    if(is.na(excludedBasins)==FALSE){
      floatData <- floatData[-which(floatData$basin %in% excludedBasins),]
    }  
    
    if(is.na(maxElevation) == FALSE){
      print("Retrieving bathymetry elevation...")
      elevation<-bathGet(floatData$latitude,floatData$longitude)
      floatData$elevation<-elevation
      floatData <- floatData %>% 
        group_by(floatID) %>% 
        mutate(coastalProfiles = sum(elevation>maxElevation,na.rm=TRUE),
               coastalPer = coastalProfiles/profiles) %>% 
        #Only keep floats where fewer than 25 percent of profiles have bathymetric
        #measurements exceeding maximum elevation allowed
        filter(coastalPer < 0.25) %>% 
        select(-coastalProfiles, -coastalPer)
    }
    
    floatData.summary<-floatData %>% 
      group_by(directory,float,floatID,year,basin,minJday,maxJday,profiles) %>% 
      summarize(n = n()) %>% 
      select(-n)
    
    write_csv(floatData.summary,"data/timeseriesData.csv",
              )
    return("Done")
  }
  
  if(timeseries==FALSE){
    print("Defining basins...")
    floatData<-floatData %>% 
      group_by(float) %>% 
      mutate(basin = basinFinder(median(longitude,na.rm=TRUE),median(latitude,na.rm=TRUE)))
    
    if(length(includedBasins)>0){
      floatData <- floatData[which(floatData$basin %in% includedBasins),]
    }
    
    if(length(excludedBasins>0)){
      floatData.test <- floatData[-which(floatData$basin %in% excludedBasins),]
    }
    
    if(is.na(maxElevation) == FALSE){
      elevation<-bathGet(floatData$latitude,floatData$longitude)
      floatData$elevation<-elevation
      
      floatData <- floatData %>% 
        filter(elevation >= maxElevation)
    }
    
    #basins<-basinFinder(floatData$longitude,floatData$latitude)
    write_csv(floatData,paste("data/profileData.csv",sep="/"),
              )
    return("Done")
  }
}

#################################################################################
#################################################################################
# getMetaData mais pour l'oxygene avec des conditions différentes puisqu'il faut également l'oxygène en Delayed Mode
# en plus de la chla et du BBP
getMetaData_DOXY<-function(params = c("CHLA","BBP700","DOXY"),
                           includedBasins = NA,
                           excludedBasins = NA, 
                           maxElevation = NA,
                           timeseries = TRUE){
  require(lubridate)
  require(stringr)
  #require(tidyverse)
  require(dplyr)
  
  #Create temporary file as tmp
  floatData <- fread("/home/master/Documents/index/argo_bio-profile_index.txt")
  
  file_split <- str_split(floatData$file,"/",simplify=T)
  floatData$directory <- file_split[,1]
  floatData$float <- file_split[,2]
  floatData$date_update <- ymd_hms(floatData$date_update)
  floatData$date<-ymd_hms(floatData$date)
  
  #Check columns for presence/absence of PRES, CHLA, BBP; add column for duration between
  #first and last profile for a given float
  for(i in 1:length(params)){
    temp<-data.frame(grepl(params[i],floatData$parameters))
    colnames(temp) = params[i]
    floatData<-cbind(floatData,temp)
    floatData<-floatData[-which(temp==FALSE),]
    file_split<-file_split[-which(temp==FALSE),]
  }
  
  floatData$year = year(floatData$date)
  floatData$jDay = yday(floatData$date)
  floatData<-floatData %>% 
    mutate(timeseries = year - min(year,na.rm = TRUE),
           floatID = paste(float,timeseries,sep=""))
  
  params_split <- str_split(floatData$parameters," ",simplify=T) #spliter parametres
  data_mode_split <- str_split(floatData$parameter_data_mode,"",simplify=T) #spliter data_mode
  
  NoCol <- vector()
  bad_DOXY <- vector()
  bad_CHLA <- vector()
  
  #selectionner les profils de CHLA qui ne sont pas en adjusted ou delayed :
  for(i in 1:nrow(floatData)){
    param_CHLA = "CHLA" %in% params_split[i,]
    if(param_CHLA) {
      NoCol[i] <- which(params_split[i,] == "CHLA")
      if (data_mode_split[i,NoCol[i]]=="D"|data_mode_split[i,NoCol[i]]=="A"){
        bad_CHLA[i] = NA
      }
      else if(data_mode_split[i,NoCol[i]]!="D"|data_mode_split[i,NoCol[i]]=="A"){
        bad_CHLA[i] <- i
      }
    }
    if(param_CHLA == FALSE) {
    }
  }
  bad_CHLA <- bad_CHLA[!is.na(bad_CHLA)]
  floatData <- floatData[-bad_CHLA,]
  data_mode_split <- data_mode_split[-bad_CHLA,]
  params_split <- params_split[-bad_CHLA,]
  
  # Il reste les profils avec CHLA = A/D et quand pas CHLA
  
  #selectionner les profils de DOXY qui ne sont pas en dalayed :
  for(i in 1:nrow(floatData)){
    param_DOXY = "DOXY" %in% params_split[i,]
    if(param_DOXY) {
      NoCol[i] <- which(params_split[i,] == "DOXY")
      if (data_mode_split[i,NoCol[i]]=="D"){
        bad_DOXY[i] = NA
      }
      else if(data_mode_split[i,NoCol[i]]!="D"){
        bad_DOXY[i] <- i
      }
    }
    if(param_DOXY == FALSE) {
      #bad_DOXY[i] <- i
    }
  }
  bad_DOXY <- bad_DOXY[!is.na(bad_DOXY)]
  floatData <- floatData[-bad_DOXY,]
  # Il reste les profils avec DOXY = D et quand pas DOXY
  
  floatData <- floatData[-which(is.na(floatData$date)),]
  
  if(timeseries==TRUE){
    print("Defining basins...")
    floatData<-floatData %>% 
      group_by(float) %>% 
      mutate(basin = basinFinder(median(longitude,na.rm=TRUE),median(latitude,na.rm=TRUE))) %>% 
      group_by(floatID) %>% 
      mutate(profiles = n(),
             minJday = min(jDay,na.rm=TRUE),
             maxJday = max(jDay,na.rm=TRUE)) %>% 
      #Only keep floats with more than 32 profiles (the minimum required to provide 
      #an annual time series with sampling resolution of 10 days, allowing for 2 missing
      #profiles at the start or end of the year)
      filter(profiles >= 32,
             minJday <= 25,
             maxJday >= 340)
    
    if(is.na(includedBasins)==FALSE){
      floatData <- floatData[which(floatData$basin %in% includedBasins),]
    }
    if(is.na(excludedBasins)==FALSE){
      floatData <- floatData[-which(floatData$basin %in% excludedBasins),]
    }  
    
    if(is.na(maxElevation) == FALSE){
      print("Retrieving bathymetry elevation...")
      elevation<-bathGet(floatData$latitude,floatData$longitude)
      floatData$elevation<-elevation
      floatData <- floatData %>% 
        group_by(floatID) %>% 
        mutate(coastalProfiles = sum(elevation>maxElevation,na.rm=TRUE),
               coastalPer = coastalProfiles/profiles) %>% 
        #Only keep floats where fewer than 25 percent of profiles have bathymetric
        #measurements exceeding maximum elevation allowed
        filter(coastalPer < 0.25) %>% 
        select(-coastalProfiles, -coastalPer)
    }
    
    floatData.summary<-floatData %>% 
      group_by(directory,float,floatID,year,basin,minJday,maxJday,profiles) %>% 
      summarize(n = n()) %>% 
      select(-n)
    
    write_csv(floatData.summary,"data/timeseriesData.csv",
              )
    return("Done")
  }
  
  if(timeseries==FALSE){
    print("Defining basins...")
    floatData<-floatData %>% 
      group_by(float) %>% 
      mutate(basin = basinFinder(median(longitude,na.rm=TRUE),median(latitude,na.rm=TRUE)))
    
    if(length(includedBasins)>0){
      floatData <- floatData[which(floatData$basin %in% includedBasins),]
    }
    
    if(length(excludedBasins>0)){
      floatData.test <- floatData[-which(floatData$basin %in% excludedBasins),]
    }
    
    if(is.na(maxElevation) == FALSE){
      elevation<-bathGet(floatData$latitude,floatData$longitude)
      floatData$elevation<-elevation
      
      floatData <- floatData %>% 
        filter(elevation >= maxElevation)
    }
    
    #basins<-basinFinder(floatData$longitude,floatData$latitude)
    write_csv(floatData,paste("data/profileData.csv",sep="/"))
    return("Done")
  }
}
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
ncImport<-function(floatList = "data/timeseriesData.csv",
                   Path = "/home/master/Documents/data/",
                   outputDirectory = "data/rawCSV/",
                   overwrite=TRUE){ 
  require(RNetCDF)
  require(zoo)
  require(curl)
  require(stringr)
  #require(tidyverse)
  require(lubridate)
  #ranger les flotteur, profils et indiquer ou ils sont (path)
  floatData <- read_csv(floatList)%>% 
    group_by(directory,float,basin) %>% 
    summarize(n=n())

  for(i in 1:nrow(floatData)){
    print(paste("File",i,"of",nrow(floatData)))
    currentFloat <- floatData$float[i]
    currentDirectory <- floatData$directory[i]
    currentBasin <- floatData$basin[i]
    
    if(overwrite==FALSE){
      existingFiles <- list.files(outputDirectory)
      if(length(grep(currentFloat,existingFiles))!=0){
        next
      }
    }
    
    floatPath = paste(Path,currentDirectory,"/",currentFloat,"/",currentFloat,"_Sprof.nc",sep="")
    
    argo <- read.nc(open.nc(floatPath)) #lire le fichier nc
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
    
    for(j in 1:length(names(argoSub))){
      outputMatrix[,j] = vectorize(argoSub[[j]],colNames[j],nCols,nRows)
    }
    
    outputDF <- data.frame(outputMatrix,stringsAsFactors = FALSE)
    colnames(outputDF)=colNames
    
    outputDF<-cbind(basin=currentBasin,float=currentFloat,cycleNumber,longitude,latitude,
                    date = ddays(jDate)+ymd_hms(refDate),direction=rep(direction,each=nRows),
                    year = year(ddays(jDate)+ymd_hms(refDate)),outputDF)
    
    outputDF<-outputDF %>% 
      filter(direction=="A",
             as.numeric(PRES)<=1000) %>% 
      mutate(floatID = paste(float,year-min(year,na.rm=TRUE)+1,sep="_")) %>% 
      select(-direction)
    
    print("Writing CSV")
    write_csv(outputDF, paste(outputDirectory,currentFloat,".csv",sep=""))
  }
}

################################################################################
################################################################################
################################################################################
#csvQC function takes a raw CSV file (generated by ncImport function) as cleans
#parameter values based on corresponding QC values, shifts dates for Southern Hemisphere
#floats, and removes incomplete profiles introduced by QC step. If unlink is set
#to TRUE, source .CSV files will be deleted after saving adjusted .CSV files.
################################################################################
################################################################################
################################################################################
#grep("2902196",list.files(inputDirectory))
csvQC <- function (inputDirectory = "data/Chla_BBP_data/rawCSV/",
                   outputDirectory = "data/Chla_BBP_data/qcCSV/",
                   overwrite = FALSE,
                   unlink = FALSE){
  library(dplyr)
  library(lubridate)
  # library(RNetCDF)
  library(stringr)
  library(readxl)

  files <- list.files(inputDirectory,pattern=".csv",full.names = T)
  noise_floats <- data.frame(noise = c(rep(NA,length(files))), float = c(rep(NA,length(files))))
  
  for(i in 1:length(files)){
    print(paste("File",i,"of",length(files)))
    import <- read_csv(files[i], show_col_types = F) %>%
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
    
    ############################################################################ 
# determination du bruit pour chaque flotteur    
    NOISE <- import %>% 
      select("PRES","date","cycleNumber","CHLA") %>%
      filter(is.na(CHLA)==FALSE) %>%
      group_by(cycleNumber) %>%
      mutate(nRows = length(cycleNumber)) %>%
      filter(nRows>11) %>%
      mutate(filtered = rollapply(CHLA, 11, min, fill = "extend")) %>%
      mutate(filtered = rollmax(filtered, 11, fill = "extend")) %>%
      filter(PRES>300) %>%
      ungroup() %>%
      mutate(bbl_noise = CHLA - filtered) %>%
      mutate(mediane = median(bbl_noise)) %>%
      
      # filter((bbl_noise<0)==F) %>%
      mutate (bin = NA)
    
    if(length(NOISE$CHLA)!=0){ # ne calcule pas de valeur si as de données

    # ranger par tranche de 50m sous 300m
    for(j in 1:length(NOISE$CHLA)){
      for(l in 1:14){
        if(NOISE$PRES[j]<(300+l*50) & NOISE$PRES[j]>=(300+(l-1)*50)){
          NOISE$bin[j] <- l
        }
      }
    }
    
    # 
    noise <- NOISE %>%
      group_by(cycleNumber, bin) %>%
      summarise(bbl_noise = bbl_noise, mediane = mediane, moyenne = mean(bbl_noise), test = moyenne < (2*mediane)) %>%
      filter(moyenne < (2*mediane)) %>% # retirer si moyenne superieure a 2 fois la mediane
      ungroup() %>%
      summarize(mediane_final = median(bbl_noise)) %>%
      as.numeric()
    
    noise_floats[i,1] <- noise
    noise_floats$float[i] <- currentFloat
    }
   ############################################################################

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
        import[pColNo][import[qcColNo]==4 | is.na(import[qcColNo])==TRUE]=NA #Set bad data (4+NA) to NA pour le realtime
        import[adjPColNo][import[adjQcColNo] == 3 | import[adjQcColNo] == 4]=NA #Set bad (3+4) adjData to NA pour le adjusted
        if(all(is.na(import[adjPColNo]))==FALSE){ #s'il y a des paramètres BBP ajustés : les déplacer dans la colonne BBP
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
        import[adjPColNo][import[adjQcColNo]==4| is.na(import[adjQcColNo])==TRUE |import[adjQcColNo]==3]=NA #Set adjParam cols with ADJ_PARAM_QC of 3 and 4 or NA to NA
        import[pColNo] = apply(import[,c(pColNo,adjPColNo)],1,mean,na.rm=TRUE)
        rm(pColNo,qcColNo,adjPColNo,adjQcColNo)
      }
      else{
        next}
    }

  import=import[-grep("_QC$|_ADJUSTED$|_ERROR$",cols)] #Remove QC, ADJUSTED and ERROR columns

  ### %m+% ensures that number of days in a month aren't exceeded when adding 6 months (as in
  ### 2004-01-31 + months(1))

  ###Create column of shifted dates; dates for Southern Hemisphere profiles will be shifted 6 months;
  ###dates for Northern Hemisphere profiles are not be changed
  ###This only changes dates for floats that only occupy the Southern Hemisphere. Floats crossing the
  ###equator will not be affected
  import$shift = NA
  import$jDay = yday(import$date)
  if(all(import$latitude<0)){
    import$date = ymd_hms(import$date)
    import$date[which(import$latitude<0)]=import$date %m+% months(6)
    import$shift = 6
  }
  if(all(import$latitude>0)){
    import$shift = 0
  }

  condensed <- import %>%
    group_by(cycleNumber) %>%
    #There are 1 or 2 floats with completely missing PSAL/TEMP data; the filter
    #settings remove all corresponding data
    filter(all(is.na(date))==FALSE,
           all(is.na(PRES))==FALSE,
           all(is.na(PSAL))==FALSE,
           all(is.na(TEMP))==FALSE,
           all(is.na(CHLA))==FALSE,
           all(is.na(BBP700))==FALSE)%>%
    ungroup()

  #enlever les profils avec moins de 5 lignes
  if(nrow(condensed)<5){
    print("Reject!")
    next}

  final <- condensed %>%
    filter(is.na(PRES)==FALSE) %>%
    dplyr::group_by(cycleNumber) %>%
    dplyr::mutate(minPres = min(PRES,na.rm=TRUE),
           maxPres = max(PRES,na.rm=TRUE),
           depths = length(which(PRES < 300)))%>%
    filter(depths>15,
           minPres < 5,
           maxPres > 300)%>%
    ungroup()


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
  
  write_csv(noise_floats, "data/Chla_BBP_data/qcCSV/recap_noise.csv")
}

################################################################################
################################################################################
################################################################################
#csvFilter function takes a list of floats as generated by the csvQC function,
#
# fonction qui permet de lire le fichiers qui recapiyule les profils de BBP avec du zoo, 
# qui ouvre les fichiers correspondants et qui filtre le zoo (filtre de Hampel avec 73m de fenetre sous 180m)

rm_zoo <- function(zoo_prof,
                   inputDirectory = "data/Chla_BBP_data/qcCSV/",
                   outputDirectory = "data/Chla_BBP_data/qcCSV/"){
  
  # avoir le fichier, le WMO, le cycle du zoo
  files <- list.files(inputDirectory,pattern=".csv",full.names = F)
  files <- str_remove(files,".csv")
  WMO <- as.character(zoo_prof$V2)
  cycle <- str_extract(zoo_prof$V4, "[1-9][0-9]*(?!.*[0-9]{3})")
  
  # pour chaque cycle, lire le fichier et appliquer la fonction get_spikes
  for(i in 1:length(cycle)){
    if(WMO[i] == "6901647" & cycle[i] == "43"){
        file <- read_csv(paste(inputDirectory,WMO[i],".csv",sep=""),show_col_types = FALSE)
        print(paste("File",i,"of",length(cycle),":",WMO[i]))
        print(paste("Cycle",cycle[i]))
        profil <- file[which(file$cycleNumber==cycle[i]),]
        x <- profil$BBP700[which(profil$PRES>180)]
        p <- profil$PRES[which(profil$PRES>180)]
        p_spikes <- p[86:160]
        x_spikes <- x[86:160]
        if(!is.character(p_spikes)){
            file <- file[-which(file$cycleNumber==cycle[i] & file$PRES %in% p_spikes & file$BBP700 %in% x_spikes),]
            write_csv(file,paste(outputDirectory,WMO[i],".csv",sep=""))
        }else{
            print(p_spikes)
        }
        next}
    if(WMO[i] %in% files){
        file <- read_csv(paste(inputDirectory,WMO[i],".csv",sep=""),show_col_types = FALSE)
        print(paste("File",i,"of",length(cycle),":",WMO[i]))
        if(cycle[i] %in% file$cycleNumber){
            print(paste("Cycle",cycle[i]))
            profil <- file[which(file$cycleNumber==cycle[i]),]
            x <- profil$BBP700[which(profil$PRES>180)]
            p <- profil$PRES[which(profil$PRES>180)]
            p_spikes <- get_spikes (x,p)$p_spikes
            
            x_spikes <- get_spikes (x,p)$x_spikes
            if(!is.character(p_spikes)){
                file <- file[-which(file$cycleNumber==cycle[i] & file$PRES %in% p_spikes & file$BBP700 %in% x_spikes),]
                write_csv(file,paste(outputDirectory,WMO[i],".csv",sep=""))
            }else{
                print(p_spikes)
            }
      }else{
        print("No corresponding cycle")
      }
    }
  }
}
    


################################################################################
################################################################################
csvAdjust <- function(inputDirectory = "data/Chla_BBP_data/qcCSV/",
                      outputDirectory = "data/Chla_BBP_data/finalCSV/",
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
        dplyr::group_by(floatID,date,basin) %>% 
        summarize_all(mean) %>% 
        select(-PRES)
      
      #Merges above vectors into a single data frame; removes missing parameter values
      #and bins by depth
      originalData <- data.frame(date = dateVector, year = yearVector, PRES = presVector,
                                 param = paramNameVector,
                                 value = paramVector, dPres = dPresVector,row.names = NULL) %>% 
        mutate(timeSeries = year - min(year)+1,
               PRES = PRES+dPres) %>% 
        dplyr::group_by(date,param) %>% 
        #Removes all NA values; NOTE! If there are no measurements for any param at 
        #less than 5 meters, this can create profiles with minimum depths of more than 5m
        #this can conflict with parameter estimates in the csvFinal script
        filter(is.na(value)==FALSE) %>% 
        dplyr::mutate(PRES = ceiling(PRES),
               n = sum(PRES<=300),
               unique = length(unique(value[PRES<=300])))%>% 
        #Removes any profile with fewer than 15 measurements remaining after removing NA
        filter(n > 15) %>% 
        ungroup() 
      
      originalData <- originalData %>%
        dplyr::group_by(date) %>% 
        #Recalculating minPres removes profiles lacking measurements at less than 5m;
        #this is redundant, given that these profiles should have been removed in 
        #bulkCsvClean, but..
        dplyr::mutate(minPres = min(PRES,na.rm=TRUE)) %>% 
        filter(minPres <= 5) %>% 
        select(-minPres) %>% 
        ungroup()
      
      originalData <- originalData %>%
        dplyr::group_by(timeSeries,date,PRES,param) %>% 
        #Averages by depth
        dplyr::summarize(value = mean(value,na.rm=TRUE)) %>% 
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
          dplyr::group_by(date,param) %>% 
          filter(length(date)>5) %>%
          dplyr::mutate(value = rollmedian(value,5,fill = "extend"),
                 value = rollmean(value,5,fill="extend"))
        
        bbpVector <- adjData %>% 
          filter(param == "BBP700") %>% 
          dplyr::group_by(timeSeries) %>% 
          dplyr::mutate(bbpThresh = quantile(value,0.995)) %>% 
          filter(value < bbpThresh) %>% 
          dplyr::group_by(date,param) %>% 
          dplyr::mutate(value = rollmedian(value,7,fill = "extend"),
                 value = rollmean(value,5,fill = "extend")) %>% 
          select(-bbpThresh)
        
        adjData <- rbind(bbpVector,nonBBP) %>% 
          spread(param,value) %>% 
          ungroup()
      }
      #If BBP is not amongst measured parameters, just spreads the existing dataVector  
      if("BBP700" %in% colnames(import)==FALSE){
        adjData<-adjData %>% 
          dplyr::group_by(date,param) %>% 
          dplyr::mutate(value = rollmedian(value,3,fill = "extend"),
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
        dplyr::group_by(date) %>% 
        #Following lines fill in any empty cells created by merge
        dplyr::mutate(float = na.fill(float,fill="extend"),
               year = na.fill(year,fill="extend"),
               jDay = na.fill(jDay,fill="extend"),
               latitude = na.fill(latitude,fill="extend"),
               longitude = na.fill(longitude,fill="extend"),
               SIG_ORIG = gsw_sigma0(PSAL_ORIG,TEMP_ORIG),
               SIG = gsw_sigma0(PSAL,TEMP)) %>% 
        dplyr::group_by(float) %>% 
        dplyr::mutate(floatID = paste(float,"_", timeSeries,sep=""))
      
      adjImport<-data.frame(adjImport)
      rm(adjOrig,originalData,adjData,import)
      adjImport$jDay = yday(adjImport$date)
      #if("DOXY" %in% colnames(adjImport)){
      # if(!all(is.na(adjImport$DOXY))){
      # adjImport$AOU = oxySat(adjImport$DOXY,adjImport$TEMP,adjImport$PSAL)
      # }else{
      #   next}#AOU
      write_csv(adjImport,
                paste(outputDirectory,currentFloat,".csv",sep=""))
      if(unlink==TRUE){
        unlink(files[i])
      }
    }
}
##############################################################################

csvAdjust_DOXY <- function(inputDirectory = "data/qcCSV/",
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
             #jDay = na.fill(jDay,fill="extend"),
             latitude = na.fill(latitude,fill="extend"),
             longitude = na.fill(longitude,fill="extend")) %>% 
             #SIG_ORIG = gsw_sigma0(PSAL_ORIG,TEMP_ORIG),
             #SIG = gsw_sigma0(PSAL,TEMP)
      group_by(float) %>% 
      mutate(floatID = paste(float,"_", timeSeries,sep=""))
    
    adjImport<-data.frame(adjImport)
    rm(adjOrig,originalData,adjData,import)
    adjImport$jDay = yday(adjImport$date)
    #if("DOXY" %in% colnames(adjImport)){
    #adjImport$AOU = oxySat(adjImport$DOXY,adjImport$TEMP,adjImport$PSAL)} #AOU
    write_csv(adjImport,
              paste(outputDirectory,currentFloat,".csv",sep=""))
    if(unlink==TRUE){
      unlink(files[i])
    }
  }
}
###############################################################################

# audit <- fread("BBP700audit112021_MF.txt", header = TRUE, sep=",")
# exceptions <- read.table("BBP700auditcoriolis012022EXCEPTIONS_MF.txt", header = TRUE, sep=",")
# colnames(audit) <- c("WMO","cycle")
# colnames(exceptions) <- c("WMO","cycle")
# 
# for(i in 1:length(exceptions$WMO)){
#   WMO <- exceptions$WMO[i]
#   cycle <- exceptions$cycle[i]
#   print(paste(WMO %in% audit$WMO,i))
#   if(WMO %in% audit$WMO){
#     if(cycle %in% audit$cycle[which(audit$WMO==WMO)]){
#       audit <- audit[-which(audit$WMO==WMO & audit$cycle==cycle),]
#     }
#   }
# }
# write_csv(audit,"new_audit.txt")


# enlever les profils mauvais de l'audit BBP
clean <- function(outputDirectory="data/finalCSV/",
                  audit=read.table("new_audit.txt", header = TRUE, sep=",")){
  
    floats <- list.files("data/Chla_BBP_data/finalCSV")
    floats <- str_replace(floats, ".csv","")
    
    WMO <- audit$WMO
    cycle <- audit$cycle
    
    for (i in 1:length(WMO)){
      if (WMO[i] %in% floats){
        print(paste(WMO[i],": File",i,"of",length(WMO),sep=" "))
        import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/",WMO[i],".csv"),show_col_types = F)
        if (cycle[i] %in% import$cycleNumber){
          print(paste("cycle",cycle[i],sep=" "))
          import <- import[-which(import$cycleNumber == cycle[i]),]
          write_csv(import, paste(outputDirectory,WMO[i],".csv",sep=""))
        }
      }else{
        next
      }
    }
    
    # enlever mauvais cycles supp
    if ("6901151" %in% floats){
      file.remove("data/Chla_BBP_data/finalCSV/6901151.csv")
    }
    
    if ("7900591" %in% floats){
      file.remove("data/Chla_BBP_data/finalCSV/7900591.csv")
    }
    
    if ("6901511" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901511.csv"))
      import <- import[-which(import$cycleNumber>=79),]
      write_csv(import, paste(outputDirectory,"6901511.csv",sep=""))
    }
    
    if ("2902092" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/2902092.csv"))
      import <- import[-which(import$cycleNumber>=121),]
      write_csv(import, paste(outputDirectory,"2902092.csv",sep=""))
    }
    
    if ("2902242" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/2902242.csv"))
      import <- import[-which(import$cycleNumber==124),]
      write_csv(import, paste(outputDirectory,"2902242.csv",sep=""))
    }
    
    if ("3901496" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/3901496.csv"))
      import <- import[-which(import$cycleNumber==160),]
      write_csv(import, paste(outputDirectory,"3901496.csv",sep=""))
    }
    
    if ("3901531" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/3901531.csv"))
      import <- import[-which(import$cycleNumber==92),]
      write_csv(import, paste(outputDirectory,"3901531.csv",sep=""),
                )
    }
    
    if ("3903586" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/3903586.csv"))
      import <- import[-which(import$cycleNumber==156),]
      write_csv(import, paste(outputDirectory,"3903586.csv",sep=""))
    }
    
    if ("5903593" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5903593.csv"))
      import <- import[-which(import$cycleNumber==43),]
      write_csv(import, paste(outputDirectory,"5903593.csv",sep=""))
    }
    
    if ("5904021" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904021.csv"))
      import <- import[-which(import$cycleNumber==161),]
      write_csv(import, paste(outputDirectory,"5904021.csv",sep=""))
    }
    
    if ("5904470" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904470.csv"))
      import <- import[-which(import$cycleNumber==47),]
      write_csv(import, paste(outputDirectory,"5904470.csv",sep=""))
    }
    
    if ("5904474" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904474.csv"))
      import <- import[-which(import$cycleNumber==146),]
      write_csv(import, paste(outputDirectory,"5904474.csv",sep=""))
    }
    
    if ("5904856" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904856.csv"))
      import <- import[-which(import$cycleNumber==33),]
      write_csv(import, paste(outputDirectory,"5904856.csv",sep=""))
    }
    
    if ("5905140" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5905140.csv"))
      import <- import[-which(import$cycleNumber==55),]
      write_csv(import, paste(outputDirectory,"5905140.csv",sep=""))
    }
    
    if ("5905982" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5905982.csv"))
      import <- import[-which(import$cycleNumber==13),]
      write_csv(import, paste(outputDirectory,"5905982.csv",sep=""))
    }
    
    if ("5906040" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5906040.csv"))
      import <- import[-which(import$cycleNumber==174),]
      write_csv(import, paste(outputDirectory,"5906040.csv",sep=""))
    }
    
    if ("6901437" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901437.csv"))
      import <- import[-which(import$cycleNumber==111),]
      write_csv(import, paste(outputDirectory,"6901437.csv",sep=""))
    }
    
    if ("6901492" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901492.csv"))
      import <- import[-which(import$cycleNumber==120),]
      write_csv(import, paste(outputDirectory,"6901492.csv",sep=""))
    }
    
    if ("6901576" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901576.csv"))
      import <- import[-which(import$cycleNumber==85),]
      write_csv(import, paste(outputDirectory,"6901576.csv",sep=""))
    }
    
    if ("6901861" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901861.csv"))
      import <- import[-which(import$cycleNumber==24|import$cycleNumber==30|import$cycleNumber==31|import$cycleNumber==37),]
      write_csv(import, paste(outputDirectory,"6901861.csv",sep=""))
    }
    
    if ("6902904" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6902904.csv"))
      import <- import[-which(import$cycleNumber==27),]
      write_csv(import, paste(outputDirectory,"6902904.csv",sep=""))
    }
    
    if ("6903026" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6903026.csv"))
      import <- import[-which(import$cycleNumber==71),]
      write_csv(import, paste(outputDirectory,"6903026.csv",sep=""))
    }
    
    if ("5904982" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904982.csv"))
      import <- import[-which(import$cycleNumber==48),]
      write_csv(import, paste(outputDirectory,"5904982.csv",sep=""))
    }

    if ("3901531" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/3901531.csv"))
      import <- import[-which(import$cycleNumber>=125 & import$cycleNumber<=139),]
      write_csv(import, paste(outputDirectory,"3901531.csv",sep=""))
    }

    if ("2902158" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/2902158.csv"))
      import <- import[-which(import$cycleNumber>=44 & import$cycleNumber<=240),]
      write_csv(import, paste(outputDirectory,"2902158.csv",sep=""))
    }
    
    if ("3901498" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/3901498.csv"))
      import <- import[-which(import$cycleNumber==285),]
      write_csv(import, paste(outputDirectory,"3901498.csv",sep=""))
    }
    
    if ("5904679" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904679.csv"))
      import <- import[-which(import$cycleNumber==23),]
      write_csv(import, paste(outputDirectory,"5904679.csv",sep=""))
    }
    
    if ("5904685" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904685.csv"))
      import <- import[-which(import$cycleNumber==84),]
      write_csv(import, paste(outputDirectory,"5904685.csv",sep=""))
    }
    
    if ("5904766" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904766.csv"))
      import <- import[-which(import$cycleNumber==108),]
      write_csv(import, paste(outputDirectory,"5904766.csv",sep=""))
    }
    
    if ("6901764" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901764.csv"))
      import <- import[-which(import$cycleNumber==207),]
      write_csv(import, paste(outputDirectory,"6901764.csv",sep=""))
    }

    if ("3901497" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/3901497.csv"))
      import <- import[-which(import$cycleNumber==37),]
      write_csv(import, paste(outputDirectory,"3901497.csv",sep=""))
    }
    
    if ("5904842" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904842.csv"))
      import <- import[-which(import$cycleNumber==88),]
      write_csv(import, paste(outputDirectory,"5904842.csv",sep=""))
    }
    
    if ("5904846" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5904846.csv"))
      import <- import[-which(import$cycleNumber==117),]
      write_csv(import, paste(outputDirectory,"5904846.csv",sep=""))
    }
    if ("5906035" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/5906035.csv"))
      import <- import[-which(import$cycleNumber==31),]
      write_csv(import, paste(outputDirectory,"5906035.csv",sep=""))
    }
    
    if ("6901174" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901174.csv"))
      import <- import[-which(import$cycleNumber==72),]
      write_csv(import, paste(outputDirectory,"6901174.csv",sep=""))
    }
    
    if ("6901764" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6901764.csv"))
      import <- import[-which(import$cycleNumber==220),]
      write_csv(import, paste(outputDirectory,"6901764.csv",sep=""))
    }
    
    if ("6903024" %in% floats){
      import <- read_csv(paste0("data/Chla_BBP_data/finalCSV/6903024.csv"))
      import <- import[-which(import$cycleNumber==116),]
      write_csv(import, paste(outputDirectory,"6903024.csv",sep=""))
    }
    
}



###############################################################################
###############################################################################

get_spikes <- function(x, p, xerr=NA, pr=NA, max_iter=NA){
  require(pracma)
  require(arules)
  # %SPIKE_ANALYSIS a point is considered as a spike if the difference between x and a 15 points Hampel
  # %       filter of x is greater than 3 times the scaled Median Absolute Deviation
  # %       also return negative spikes separately
  # %
  # % INPUTS (units do not matter):
  #   %     x <Nx1 double> profile of bio-optical data, for example:
  #   %         - particulate backscattering (beta, bbp)
  # %         - fluorescent dissolved organic matter (fdom)
  # %         - chlorophyll a fluorescence  (fchl)
  # % OPTIONAL INPUT:
  #   %     p <Nx1 double> pressure of profile (add sampling frequency check)
  # %     xerr <1x1 double> absolute threshold for spike detection
  # %         default: 3x the scaled Median Absolute Deviation (MAD) of entire profile (non-parametric)
  # %         recommended: for fdom: ???; for bbp: ???; for fchl: ???;                 (parametric)
  # %     pr <1x1 double> minimum sampling resolution threshold
  # %         default: 3x the median sampling resolution of entire profile             (non-parametric)
  # %         recommended: 5 dBar (typically no spike layers are found below 2.25 dBar)(parametric)
  # %     max_iter <1x1 integer> maximum number of intereation to run Hampel Filter
  # %         set to 0 to run the Hampel filter only once
  # %         default: 100
  # %     
  # %
  # % OUTPUTS:
  #   %     spikes <NX1 boolean> return true for every indices containing a positive spike
  # %     neg_spikes <NX1 boolean> return true for every indices containing a negative spike
  # %     spikes_p <NX1 boolean> pressure of each positive spikes
  # %     neg_spikes_p <NX1 boolean> pressure of each negative spikes
  # %
  # % References:
  #   %   Rousseeuw, P. J., and C. Croux (1993), Alternatives to the Median Absolute Deviation,
  # %     Journal of the American Statistical Association, 88(424), 1273?1283.
  # %
  # % author: Nils Haentjens
  # % created: May 1, 2018
  x <- tibble(prof = p, x = x) %>% arrange(prof) %>% pull(x)
  p <- tibble(prof = p, x = x) %>% arrange(prof) %>% pull(prof)
  
  DEBUG = F
  
  # Check input
  #nargin <- length(as.list(match.call())) -1
  nargin = 2
  
  if (nargin < 2 | length(p) == 0){p_flag = F
  }else {p_flag = T}
  
  if (nargin < 3 | length(xerr) == 0){ xerr_flag = F
  xerr <- NA
  }else {xerr_flag = T}
  
  if (nargin < 4 | length(pr) == 0){ pr_flag = F
  }else {p_flag = T} 
  #ignored if ~p_flag
  
  if (nargin < 5 | length(max_iter) == 0){max_iter = 100}
  
  
  if (!is.vector(x)) {warning('x is not a row vector (Nx1).')}
  if (p_flag & any(length(p) != length(x))) {warning('Input vectors are different size.')}
  
  # % Prepare input
  spikes <- logical(length = length(x))  #attention il faut que toutes les valeurs soient = à FALSE
  
  if (p_flag==T){
    # ignore nan of x and p for processing
    sel = which(!is.na(x) & !is.na(p) & !is.infinite(x) & !is.infinite(p))
    
    
    if (length(sel) < 3) {
      warning('GET_SPIKES:EMPTY_INPUT', 'get_spikes: Not enough valid values.')
      
      if (nargout == 2 & !p_flag) {varargout = spikes}
      else {varargout <- list(vector(), spikes, vector())}
      return(varagout)
    }
    
    p1 = p[sel]
    # ignore low sampling resolution of profile
    delta <- matrix(nrow=length(p1), ncol=1)
    delta[1,1] <-  abs(diff(p1[1:2]))
    delta[2:(length(p1)-1),1] <- abs((p1[3:length(p1)] - p1[1:(length(p1)-2)]) / 2)
    delta[length(p1),1] <- abs(diff(p1[(length(p1)-1):length(p1)]))  # / dBar 
    if (!pr_flag){ 
      pr = 3 * median(delta)
      subsel = delta < pr
      sel = sel[subsel]
    }
    if (DEBUG){
      #subplot(1,3,1)
      split.screen(c(1,3)) 
      plot(delta, p)
      set(gca, 'ydir', 'reverse') #renverser l'axe Y
    }
    # % keep only relevant sampling
    x = x[sel]
    if(length(p)!=length(sel)){
      p = p[sel]
    }else{
      p=p1
    }
  }else{
    # % only x
    sel = which(!is.na(x) & !is.infinite(x))
    x = x[sel]
  }
  
  if(length(x)>=(35*2)){
    # % Apply Hampel filter (same as median filter except onlx change values of spikes)
    xd <- hampel(x, 35) 
    # % Recursive Hampel filter
    xdo <- x
    i <- 0
    while (any(xd$y != xdo) & i < max_iter){
      xdo <- xd$y
      xd <- hampel(xd$y, 35)
      i <- i + 1
    }
    if (i == max_iter){warning('GET_SPIKES:MAX_ITER', 'get_spikes: Maximum iteration (%d) reached.', max_iter)}
    
    if (!xerr_flag){
      # % Compute Scaled Median Absolute Deviation (MAD)
      smad = -1/(sqrt(2)*erfcinv(3/2)) * median(abs(x-median(x)))
      # Noisy signal (no significant spike)
      if (smad == 0){
        # % format output & return
        if (nargout == 2 & !p_flag){varargout = {spikes}}
        else {varargout = list(vector(), spikes, vector()) } 
        return(varargout)
      }
      # % Set spike detection threshold to 3 scaled MAD
      xerr = 3 * smad
    }
    
    # % Get spikes
    x_delta <-  x - xd$y
    spikes[sel] = x_delta > xerr
    
    # % Set output
    if (p_flag){
      spikes <- which(spikes[sel]==T)
      
    }
    
    if(length(spikes)!=0){
      x_spikes <- x[spikes] # profil déspiké
      p_spikes <- p[spikes] # profondeurs ou il y a les spikes
    }
    return(list(p_spikes=p_spikes, x_spikes=x_spikes ))
    
  }else{
    return(list(p_spikes="sample is too small", x_spikes="sample is too small"))
  }
}
    
###############################################################################
###############################################################################
# supprimer le profil si pas de latitude ou longitude
for (i in 1:8){
    files <- list.files("Arabian/final",pattern=".csv",full.names = T)
    print(paste("File",i,"of",length(files)))
    import <- read_csv(files[i],show_col_types = F)
    cycles <- unique(import$cycleNumber)
    currentFloat <- import$float[1]
    for (j in 1:length(cycles)){
        profil <- import[which(import$cycleNumber == cycles[j]),]
        if(all(is.na(profil$latitude)) | all(is.na(profil$longitude))){
          import <- import[-which(import$cycleNumber == cycles[j]),]
        }
    }
    write_csv(import,paste("Arabian/",currentFloat,".csv"))
}

################################################################################
# créer une colonne et convertir BBP en POC avec équation de Raph 2022

bbp_POC <- function(){
  
  files <- list.files("data/finalCSV",pattern=".csv",full.names = T)
  
  for (i in 1:length(files)){
    print(paste("File",i,"of",length(files)))
    import <- read_csv(files[i],show_col_types = F)
    import <- mutate(import, POC = 38687.27*(import$BBP700)^(0.95)) #POC in mg/m3 Sauzede, 2022
    currentFloat <- import$float[1]

    write_csv(import,paste("data/Chla_BBP_data/files_POC/",currentFloat,".csv",sep=""))
    
  }

}

bbp_POC()

################################################################################

# masquer les donnees pour ne garder que le POC associé à assez de chla pour faire partie de la biomasse phyto
# au lieu du POC utiliser Cphyto et changer le seuil de 200m pour le 90e percentile 

mask_Uchida_2019 <- function(){
  
  files <- list.files("data/Chla_BBP_data/files_POC/",pattern=".csv",full.names = T)
  
  for (i in 1:length(files)){
    
    print(paste("File",i,"of",length(files)))
    import <- read_csv(files[i],show_col_types = F)
    currentFloat <- import$float[1]
  
    #below_Ze <- filter(import, PRES>max(MLD_003,ISO_1,na.rm=T))
    
    below_Ze <- filter(import, PRES>200)
    percentile_90 <- as.numeric(quantile(below_Ze$CHLA_recal,0.9, na.rm=T))
    
    CHL_masked <- import$CHLA_recal
    Cphyto_masked <- import$Cphyto
    
    if(any(CHL_masked<percentile_90)){
      range <- which(CHL_masked<percentile_90)
      CHL_masked[range] <- 0
      Cphyto_masked[range] <- 0
    }
    
    import <- mutate(import, CHL_masked = CHL_masked)
    import$Cphyto_masked <- Cphyto_masked 
 
    write_csv(import,paste("data/Chla_BBP_data/finalCSV/",currentFloat,".csv",sep=""))
    
  }
}

mask_Uchida_2019()
    
##################################################################################

bbp_POC_Koestner <- function(){
  
  files <- list.files("data/Chla_BBP_data/finalCSV/",pattern=".csv",full.names = T)
  k1 <- 89.423
  k2 <- 0.1881
  k3 <- 0.7591
  k4 <- 0.1934
  POCmin <- 33.4
  e1 <- 1.636
  e2 <- 21.2
    
  for(i in 271:348){
    print(paste(i, "sur", length(files)))
    import <- read_csv(files[i],show_col_types = F)
    # POC Sauzède 2022
    import <- mutate(import, POC_Sauzede = 38687.27*(import$BBP700)^(0.95)) #POC in mg/m3 Sauzede, 2022
    
    # POC Koestner 2022
    zeta = import$CHLA_recal / import$BBP700
                     
    # if(zeta < 1){
    #   zeta=1
    # }
    
    POC_star <- k1 * import$BBP700^k2 * zeta^k3 * zeta^(k4*log10(import$BBP700))
    
    import <- mutate(import, POC_Koestner = POC_star, 
                     zeta = zeta)
    currentFloat <- import$float[1]
    
    # fonction de biais pour corriger quand POC<POCmin
    for (j in 1:length(import$BBP700)){
      
        if(is.nan(POC_star[j]) | is.na(POC_star[j])){
          # print(paste("a donné un NaN", currentFloat))
          import$POC_Koestner[j] <- 0
          next
        }
      
        if((POC_star[j]<POCmin)){
          # print("<POCmin")
          import$POC_Koestner[j] <- e1 * POC_star[j] - e2
          if(import$POC_Koestner[j]<0){
            # print(paste("a donné un résultat <0", currentFloat))
            import$POC_Koestner[j] <- 0
          }
        }
      
    }
    
    
    write_csv(import,paste("data/Chla_BBP_data/files_POC/",currentFloat,".csv",sep=""))
  }
  
}

bbp_POC_Koestner()

######################################################################
Cphyto <- function(){
  
  files <- list.files("data/Chla_BBP_data/finalCSV/",pattern=".csv",full.names = T)
  
  for (i in 1:length(files)){
    print(paste("File",i,"of",length(files)))
    import <- read_csv(files[i],show_col_types = F)
    BBP700 <- import$BBP700

    import <- mutate(import, BBP470 = BBP700 * (470/700)^(-0.78), # Uchida et al. 2020
                             Cphyto = 12128 * BBP470 + 0.59) # Graff et al 2015
    
    currentFloat <- import$float[1]
    
    write_csv(import,paste("data/Chla_BBP_data/files_POC/",currentFloat,".csv",sep=""))
    
  }
  
}

Cphyto()



files <- c("5903274_2","5903274_3","5903274_4","5903714_2","5903714_3","5903714_4","5904021_2","5905988_2","5906293_2")

files <- c("5903887_4","5904479_2","6901486_2","6901486_3","6901486_4","6901519_2","6901523_2","6902547_2","6903549_2","6903549_3","6903550_2")

inputMatrix <- read_csv("data/all/DOXY_data_all.csv") 

for (i in 1:5){

  inputMatrix_2 <- inputMatrix %>% 
    dplyr::filter(PARAM=="DOXY" & cluster==i) %>%
    # dplyr::select("jDay","PRES","PARAM",value = files[i]) %>%
    dplyr::select("jDay","PRES","PARAM","value") %>%
    dplyr::group_by(PRES) %>%
    dplyr::summarise(PRES=unique(PRES), mean_value = mean(value), sd = sd(value))
  
  a <- ggplot(inputMatrix_2, aes(PRES,mean_value)) + 
    geom_line() +
    geom_ribbon(aes(ymin=mean_value-sd,ymax=mean_value+sd), alpha=0.2) +
    scale_x_reverse(expand=c(0,0),limits=c(1000,0)) + 
    scale_y_continuous(expand=c(0.01,0.01)) +
    coord_flip() + 
    xlab("pression") +
    ylab("[O2] µmol/kg") +
    theme_bw() +
    ggtitle(paste("cluster",i)) +
    ylim(0,310)
  
  # b <- ggplot(inputMatrix_2, aes(PRES,mean_value)) + 
  #   geom_line() +
  #   geom_ribbon(aes(ymin=mean_value-sd,ymax=mean_value+sd), alpha=0.2) +
  #   scale_x_reverse(expand=c(0,0),limits=c(1000,0)) + 
  #   scale_y_continuous(expand=c(0.01,0.01)) +
  #   coord_flip() + 
  #   xlab("pression") +
  #   ylab("[O2] µmol/kg") +
  #   theme_bw() +
  #   ylim(0,50)
  # 
  # c=a+b
  print(a)
  
  # ggsave(plot=c,paste("figures/interpolatedTimeseries/O2/AN/profiles/",files[i],".png",sep=""), height=3, width=4, dpi=300)
  
}

############################################
recal_CHLA <- function(){

noise <- read_csv("data/Chla_BBP_data/recap_noise.csv")
files <- list.files("data/Chla_BBP_data/files_POC/",pattern=".csv",full.names = T)

for (i in 146:length(files)){
  print(paste("File",i,"of",length(files)))
  import <- read_csv(files[i],show_col_types = F)
  currentFloat <- unique(import$float)
  import <- mutate(import, CHLA_recal = CHLA)
  
  ## si valeurs negatives alors mettre la valeur min du profil à la valeur du 
  ## bruit et décaler tout le profil
  if(any(import$CHLA<0, na.rm=T)){
    cycleNumbers <- unique(import$cycleNumber)
    
    for(j in 1:length(cycleNumbers)){
      cycle <- import %>%
        filter(cycleNumber == cycleNumbers[j])
      if(any(cycle$CHLA<0, na.rm=T)){
        CHLA <- cycle$CHLA
        float_noise <- noise %>%
          filter(float == currentFloat) %>%
          select(noise) %>%
          as.numeric()
        if(is.na(float_noise)==F){
          import$CHLA_recal[which(import$cycleNumber==cycleNumbers[j])] <- CHLA + abs(min(CHLA, na.rm = T)) + abs(float_noise)
          # print(paste("cycle :",cycleNumbers[j]))
        }
        if(is.na(float_noise)==T){
          import$CHLA_recal[which(import$cycleNumber==cycleNumbers[j])] <- CHLA
          # print(paste("cycle :",cycleNumbers[j]))
        }
      }
    }
    
  }
  write_csv(import,paste("data/Chla_BBP_data/finalCSV/",currentFloat,".csv",sep=""))
}


}

recal_CHLA()

