#SET WORKING DIRECTORY TO PATH OF argoSVD directory:
#setwd("/Users/nicholasbock/Projects/argo/argoSVD/")

require(tidyverse)
require(lubridate)
require(zoo)
require(SolarSpectrum)

Sys.setenv(tz = "UTC")
options(warn = 1)


###############################################################################
###############################################################################
###################SCRIPTS SUPPORTING PARAMETER DETERMINATION##################
###############################################################################
###############################################################################
par_prof_cornec<-function(chl_profile,depth_profile,lon,lat,time_posix,month,Table_KD,table_coeff_corrlight,region) {
  #browser()
  #Calculate the PAR profiles using the Kd estimated from bio-optical relationships between the chla attenuation at 490 nm
  # Inputs : chl_profile
  #         : depth_profile
  #         : lon
  #         : lat
  #         : time_posix = time in a posixtc format 
  #         : month of the profile = numeric variable (1-12)
  #         : Table_KD = table with the coefficients a and b estimated for the global irradiance/chla BGC argo datatbase
  #         : table_coeff_corrlight = table with the cloud cover corrections estimated at each 10? lat band and per month
  #         : region = acronym of the region of the profile (from Zone.R function)
  # 
  # Output : PAR_cor : a vector with the PAR profiles estimated at regular depth bewteen 0 and 500m
  
  # Required : SolarSpectrum = model of surface irradiance from clear sky, Gregg and Carder, 1990
  
  library(SolarSpectrum)
  
  
  # time as.posixt
  
  
  
  
  # interpolate Fchla per meter
  # dep_int<-c(0:max(depth_profile,na.rm=TRUE))
  # chl_int<-NA
  dep_int<-depth_profile
  chl_int<-chl_profile
  
  
  # attribute the bio-optical coefficient according to the region of the profile
  if(region %in% c("NASTG","NPSTG","SASTG","SPSTG","ARCH")) {region<-"DAZ"}
  if(region %in% c("CHINS","ASEW","PSEW","WMS","AUS","KURIO","IOMZ","IEQ","MOONS")) {region<-"DBZ"}
  if(region %in% c("ASZ_SIZ","NAC","STZ","PFZ","SAZ")) {region<-"GOHZ"}
  if(region %in% c("NASPG","NS","NPSPG","ARCT")) {region<-"SHAZ"}
  
  
  # Kd 490 MOREL 2007: original coefficients and relation
  # KD490<-NA
  # KD490<- 0.0166 + (0.0773*(chl_int^0.6715))
  
  
  #retrieve a and b coefficient according to the region
  a_coeff<-NA
  b_coeff<-NA
  
  a_coeff<-as.numeric(paste(Table_KD$a[Table_KD$i==region]))
  b_coeff<-as.numeric(paste(Table_KD$b[Table_KD$i==region]))
  
  # if no coeff attributed, use coeffs of Morel 2007
  if(length(a_coeff)==0|length(b_coeff)==0){
    a_coeff<-0.0773
    b_coeff<-0.6715
  }
  
  #retrieve KD490 profiles
  KD490<-NA
  KD490<- 0.0166 + (a_coeff*(chl_int^b_coeff))
  
  
  #KD PAR
  KD_PAR<-NA
  KD_PAR<- 0.0665 + 0.874*KD490 - (0.00121*(KD490^-1))
  
  
  
  # KD from Kim et al., 2005, other version of Kd relationship with chla profile
  # KD_PAR<-NA     
  # KD_PAR<- 0.0232 + (0.074*(chl_int^0.674))
  
  KD_PAR[which(is.na(KD_PAR)==T)]<-0
  
  # integrate the Kd at depth
  KD_sum<-NA
  KD_sum<-cumsum(KD_PAR)
  
  # cloud cover coeff
  coeff<-NA
  coeff<-table_coeff_corrlight[month,which.min(abs(seq(-70,70,10)-lat))]
  
  if(is.infinite(coeff) | is.na(coeff) | coeff > 1){coeff=1}
  
  ##instantaneous PAR at noon (from the solarspectrum package) corrected from the cloud coverage
  par0_th<-NA
  par0_th=SolarSpectrum.PARi(lon,lat,time_posix)[[2]]*(10^6)*coeff
  
  #PAR profile
  PAR_cor<-NA
  PAR_cor<-par0_th*(exp(-KD_sum))
  
  
  return(data.frame(PRES = dep_int,
                    PAR_EST = PAR_cor))
  
}

daily_PAR_given_depth<-function(lon,lat,time_posix,par_profile,depth_profile, given_depth) {
  #browser()
  if(length(given_depth)==0){return(NA)}
  # Calculate the daily integrated PAR (over the day length) at a given depth
  
  # Inputs : par_profile (measured or estimated)
  #         : depth_profile
  #         : lon
  #         : lat
  #         : time_posix = time in a posixtc format 
  #         : given depth= the depth at which the daily PAR should be estimated
  # 
  # Output : daily_PAR_depth : the daily PAR value estimated at the given depth
  
  # Required : SolarSpectrum = model of surface irradiance from clear sky, Gregg and Carder, 1990
  
  library(SolarSpectrum)
  
  
  par0_noon<-NA
  par0_th<-NA
  par0<-NA
  Kd<-NA
  corr.factor=NA
  
  daily_par_depth<-NA
  
  #estimation of the surface instanteanous PAR at the geographical position and time from the model
  par0_th=SolarSpectrum.PARi(lon,lat,time_posix)[[2]]*10^6
  
  # estimation of the surface max PAR value from the PAR profile
  par0_noon=max(par_profile[1:5],na.rm=T)
  
  
  #calculation of a correction coefficient between the par surface and the model
  if(is.infinite(par0_noon)){par0_noon=NA}
  if (par0_th>0 & !is.na(par0_noon))
  {
    corr.factor=par0_noon/par0_th
    if(is.infinite(corr.factor) | is.na(corr.factor)){corr.factor=1}
    
    # calculation of the daily surface PAR
    par0=SolarSpectrum.PAR(lon,lat,time_posix)[[2]]*corr.factor
  }
  
  
  
  if(is.na(given_depth)==F) {
    ratio_given_depth<-NA
    # define the ratio of attenuation between the surface and the given depth from the par profile
    ratio_given_depth<-par_profile[which.min(abs(c(0:300)-given_depth))]/
      par0_noon
    # estimation of the daily PAR at the given depth
    daily_par_depth<-par0*ratio_given_depth
  }
  
  return(daily_par_depth)
}

avgMldPAR <- function(lon,lat,time_posix,par_profile,depth_profile,MLD_depth){
  #browser()
  if(min(depth_profile,na.rm=TRUE)>MLD_depth){
    next
    # depth_profile <- approx(depth_profile,par_profile,1:1000)$y
    # par_profile <-  approx(depth_profile,par_profile,1:1000)$x
  }
  depths <- c(0,depth_profile[depth_profile<=MLD_depth])
  depthDiffs <- diff(depths)
  intParProfile = c()
  for(i in 1:(length(depths)-1)){
    iPAR <- daily_PAR_given_depth(lon,lat,time_posix,par_profile,depth_profile,depth_profile[i])
    intParProfile = c(intParProfile,iPAR)
  }
  totalPAR = sum(intParProfile * depthDiffs)
  return(totalPAR/MLD_depth)
}
###############################################################################
###############################################################################
########################FUNCTIONS TO DETERMINE PARAMETERS######################
###############################################################################
###############################################################################
parEstimator <- function(inputDirectory = "data/finalCSV/",
                         outputDirectory = "data/parData/",
                         plot=TRUE){
  require(patchwork)
  require(lubridate)
  require(zoo)
  require(tidyverse)
  tableKD = read.table("table_KD.txt")[-1,]
  fileList <- list.files(inputDirectory,pattern=".csv")
  coefCor = read.table("table_coeff_corr_light.txt")
  colnames(coefCor)<-c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70)
  for(i in 1:length(fileList)){
    print(i)
    
    csvFile<-read_csv(paste(inputDirectory,fileList[i],sep=""),show_col_types = F)
    
        currentFloat <- csvFile$float[1]
        if("CHLA_ORIG" %in% colnames(csvFile)==FALSE){
          next
        }
        
        input<-csvFile %>% 
          arrange(date,PRES) %>% 
          select(floatID,basin,province,date,latitude,longitude,
                 PRES,CHLA_ORIG)
        
        if("DOWNWELLING_PAR_ORIG" %in% colnames(csvFile)){
          input<-csvFile %>% 
            arrange(date,PRES) %>% 
            select(floatID,basin,province,date,latitude,longitude,
                   PRES,DOWNWELLING_PAR_ORIG,CHLA_ORIG)
        } 
        
        outputData <- data.frame()
        for(j in 1:length(unique(input$date))){
          print(j)
          origTemp <- NA
          range <- which(input$date==unique(input$date)[j])
          currentDate <- as.POSIXct(unique(input$date)[j])
          #currentCycle <- unique(input$cycleNumber)[range][1]
          currentFloatID <- input$floatID[range][1]
          currentProvince <- input$province[range][1]
          currentBasin <- input$basin[range][1]
          currentLat <- input$latitude[range][1]
          currentLon <- input$longitude[range][1]
          
          #Interpolate DOWNWELLING_PAR_ORIG, if present; this will not fill missing starting values
          if("DOWNWELLING_PAR_ORIG" %in% colnames(input) &
             all(is.na(input$DOWNWELLING_PAR_ORIG[range]))==FALSE){
            input$DOWNWELLING_PAR_ORIG[range][input$DOWNWELLING_PAR_ORIG[range]<0]=0
            origTemp <- approx(input$PRES[range],
                               input$DOWNWELLING_PAR_ORIG[range],
                               1:max(input$PRES[range],na.rm=TRUE))$y
          }
    
          #If province is "CAL" (not represented in integrated PAR function), or if there are 
          #no CHLA measurements or province values within the range, skip profile
          if(input$province[range][1]=="CAL"|
             all(is.na(input$CHLA_ORIG[range]))==TRUE|
             all(is.na(input$province[range]))==TRUE){
            next
          }
          
          input$CHLA_ORIG[range][input$CHLA_ORIG[range]<0]=0
          
          chlTemp <- approx(input$PRES[range],
                            input$CHLA_ORIG[range],
                            1:max(input$PRES[range],na.rm=TRUE))$y
          presTemp <- approx(input$PRES[range],
                             input$CHLA_ORIG[range],
                             1:max(input$PRES[range],na.rm=TRUE))$x
          chlTemp<-na.fill(chlTemp,"extend")
          
          parEst <- par_prof_cornec(chl_profile = chlTemp,
                                    depth_profile = presTemp,
                                    lon = currentLon,
                                    lat = currentLat,
                                    time_posix = currentDate,
                                    month = month(currentDate),
                                    Table_KD = tableKD,
                                    table_coeff_corrlight = coefCor,
                                    region = as.character(currentProvince))
          
          #Set DOWNWELLING_PAR column to PAR_EST and DOWNWELLING_PAR_ORIG to NA
          profileData <- data.frame(floatID = currentFloatID,
                                    date = currentDate,
                                    #cycleNumber = currentCycle,
                                    province = currentProvince,
                                    lat = currentLat,
                                    lon = currentLon,
                                    PRES = parEst$PRES,
                                    CHLA = chlTemp,
                                    DOWNWELLING_PAR_EST = parEst$PAR_EST,
                                    DOWNWELLING_PAR_MEAS = NA,
                                    DOWNWELLING_PAR = parEst$PAR_EST,
                                    PAR_FLAG = 0)
          
          #If DOWNWELLING_PAR_ORIG is available, and if the first measurement is at 1m,
          #set DOWNWELLING_PAR to DOWNWELLING_PAR_ORIG (and PAR_FLAG to 1)
          if("DOWNWELLING_PAR_ORIG" %in% colnames(input) &
             is.na(origTemp[1])==FALSE &
             (any(origTemp[1:5]>origTemp[1])==FALSE)){
            profileData$DOWNWELLING_PAR_MEAS=origTemp
            profileData$DOWNWELLING_PAR = origTemp
            profileData$PAR_FLAG = 1
          }
          
          outputData <- rbind(outputData,profileData)
     
          if(plot == TRUE){
            output = ggplot()+
              geom_point(aes(parEst[,1],parEst[,2])) +
              theme_bw() +
              coord_flip() +
              scale_x_reverse(limits = c(300,0)) + 
              labs(x = "Depth (m)", y = "PAR")
            
            if("DOWNWELLING_PAR_ORIG" %in% colnames(input) &
               all(is.na(input$DOWNWELLING_PAR_ORIG[range]))==FALSE){
              
              output = output +
                geom_point(aes(presTemp,origTemp),color = "red", alpha = 0.3) 
            }
            
            ggsave(plot = output,
                   paste("figures/parProfiles_Est+Meas/",
                         currentFloat,"_",i,"_",j,".png",sep=""))
          }
          }
        
        #CAL floats will skip all profiles, returning an empty data frame that will 
        #subsequently create errors when iterating through files. To avoid this, don't
        #save any data frames with 0 rows.
        if(nrow(outputData)==0){
          next
        }
        
        if(all(is.na(outputData$DOWNWELLING_PAR_MEAS))==FALSE &
           plot == TRUE){
          outputPlot <- ggplot(outputData) + 
            geom_point(aes(DOWNWELLING_PAR_EST,DOWNWELLING_PAR_MEAS)) + 
            geom_abline(aes(slope = 1, intercept = 0), col = "red",size = 1,linetype = 2) + 
            theme_bw() + 
            facet_wrap(~month(date),scales="free") + 
            theme(plot.title=element_text(hjust = 0.5)) +
            ggtitle(paste(currentProvince,"|||",currentFloat))
          
          ggsave(paste("figures/parameterTestPlots/estParMeasPar/",
                       currentProvince,"_",currentFloat,".png",sep=""),
                 dpi = 300)
        }
        
        write.csv(outputData,
                  file=paste(outputDirectory,currentFloat,".csv",sep=""),
                  row.names=FALSE)
        rm(outputData)
  }
}

intPAR <- function(fileList,parFiles,plot=TRUE,
                   inputDirectory = "data/finalCSV/",
                   outputDirectory = "data/finalCSV/"){
  library(lubridate)
  library(tidyverse)
  
  tableKD = read.table("supportingFiles/table_KD.txt")[-1,]
  coefCor = read.table("supportingFiles/table_coeff_corr_light.txt")
  colnames(coefCor)<-c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70)
  for(i in 1:length(fileList)){
    print(paste(i,fileList[i]))
    
    temp<-read.csv(paste(inputDirectory,fileList[i],sep=""))
    currentFloat = temp$float[1]
  
    temp$NCLINE_PAR = NA
    temp$CAN_NCLINE_PAR = NA
    temp$CAN_PCLINE_PAR = NA
    temp$CAN_SCLINE_PAR = NA
    temp$MLD_MED_PAR = NA
    #temp$DCM_PAR = NA
    
    if(length(grep(currentFloat,parFiles))==0|
       #length(unique(temp$date))<20|
       "CHLA"%in% colnames(temp)==FALSE|
       "BBP700"%in% colnames(temp)==FALSE){
      print("insufficient data")
      write.csv(temp, file = fileList[i],
                row.names=FALSE)
      next
    }
    
    parData<-read.csv(parFiles[grep(currentFloat,parFiles)]) %>% 
        arrange(floatID,date,PRES)
    
    maxPar = max(parData$DOWNWELLING_PAR,na.rm=TRUE)
    
    for(j in 1:length(unique(temp$date))){
      print(j)

      range = which(temp$date==unique(temp$date)[j])
      parRange = which(as.character(parData$date)==as.character(unique(temp$date)[j]))
      #Set parTemp as estimated par values within parRange (values already interpolated to 1m)
      
      #Set presTemp as estimated PRES values within parRange (values already interpolated to 1m)
      presTemp = parData$PRES[parRange]
      parTemp = parData$DOWNWELLING_PAR[parRange]

      currentDate <- as.POSIXct(unique(temp$date)[j])
      currentBasin <- temp$basin[range][1]
      currentMLD <- temp$MLD_003[range][1]
      currentProvince <- temp$province[range][1]
      
      if(all(is.na(parTemp))|
         all(is.na(temp$province[range]))){
        next
      }
      
      tempPARprof <- approx(presTemp,parTemp,1:max(presTemp))$y
      tempPRESprof <-  approx(presTemp,parTemp,1:max(presTemp))$x
      
      mldMedian <- floor(median(tempPRESprof[tempPRESprof <= currentMLD],na.rm=TRUE)) 
      
      temp$MLD_MED_PAR[range] = daily_PAR_given_depth(lon = temp$longitude[range][1],
                                                      lat = temp$latitude[range][1],
                                                      time_posix = currentDate,
                                                      par_profile = tempPARprof,
                                                      depth_profile = tempPRESprof,
                                                      given_depth = mldMedian)
      
      temp$NCLINE_PAR[range] = daily_PAR_given_depth(lon = temp$longitude[range][1],
                                                     lat = temp$latitude[range][1],
                                                     time_posix = currentDate,
                                                     par_profile = tempPARprof,
                                                     depth_profile = tempPRESprof,
                                                     given_depth = temp$NCLINE_DEPTH[range][1])
      
      temp$CAN_NCLINE_PAR[range] = daily_PAR_given_depth(lon = temp$longitude[range][1],
                                                         lat = temp$latitude[range][1],
                                                         time_posix = currentDate,
                                                         par_profile = tempPARprof,
                                                         depth_profile = tempPRESprof,
                                                         given_depth = temp$CAN_NCLINE_DEPTH[range][1])
      
      temp$CAN_PCLINE_PAR[range] = daily_PAR_given_depth(lon = temp$longitude[range][1],
                                                         lat = temp$latitude[range][1],
                                                         time_posix = currentDate,
                                                         par_profile = tempPARprof,
                                                         depth_profile = tempPRESprof,
                                                         given_depth = temp$CAN_PCLINE_DEPTH[range][1])
      
      temp$CAN_SCLINE_PAR[range] = daily_PAR_given_depth(lon = temp$longitude[range][1],
                                                         lat = temp$latitude[range][1],
                                                         time_posix = currentDate,
                                                         par_profile = tempPARprof,
                                                         depth_profile = tempPRESprof,
                                                         given_depth = temp$CAN_SCLINE_DEPTH[range][1])
      
    # if(plot==TRUE){
    #   output<-ggplot() +
    #     geom_point(aes(tempPRESprof,tempPARprof)) +
    #     scale_x_reverse(limits=c(500,0)) +
    #     scale_y_continuous(limits = c(0,maxPar))+
    #     geom_vline(aes(xintercept = temp$CAN_NCLINE_DEPTH[range][1]),col="red") +
    #     geom_vline(aes(xintercept = mldMedian)) +
    #     geom_text(aes(x = (mldMedian + 10),
    #                   y = (maxPar - 20),
    #                   label = round(temp$MLD_MED_PAR[range][1],2)),
    #               col = "black")+
    #     geom_text(aes(x = (temp$CAN_NCLINE_DEPTH[range][1] + 10),
    #                   y =(maxPar - 20),
    #                   label = round(temp$CAN_NCLINE_PAR[range][1],2)),
    #               col = "red")+
    #     coord_flip() +
    #     theme_bw() +
    #     ggtitle(currentDate)
    #   ggsave(paste(parameterTestPlots/intPAR/",
    #                currentBasin,"_",currentFloat,"_",j,".png",sep=""),
    #          plot = output)
    # }

    }
    if(plot==TRUE){
      summary <- temp %>% 
        select(jDay,MLD_MED_PAR,CAN_NCLINE_PAR) %>% 
        group_by(jDay) %>% 
        summarize_all(mean,na.rm=TRUE)
      
      output<-ggplot(summary) + 
        geom_point(aes(jDay,(MLD_MED_PAR))) + 
        geom_point(aes(jDay,(CAN_NCLINE_PAR)),col="red") + 
        theme_bw() 
      
      ggsave(paste("figures/parameterTestPlots/intPAR_TS/",
                   currentBasin,"_",currentFloat,".png",sep=""),
             plot = output)
    }
    
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }
}

intMldPAR <- function(fileList, replace = FALSE,
                      inputDirectory = "data/finalCSV/",
                      outputDirectory = "data/finalCSV/"){
  library(lubridate)
  library(SolarSpectrum)
  tableKD = read.table("supportingFiles/table_KD.txt")
  coefCor = read.table("supportingFiles/table_coeff_corr_light.txt")
  colnames(coefCor)<-c(-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70)
  for(i in 1:length(fileList)){
    print(paste(i,fileList[i]))
    
    temp<-read.csv(paste(inputDirectory,fileList[i],sep=""))
    
    temp$MLD_INT_PAR = NA
    
    if(length(unique(temp$date))<20|
       "CHLA"%in% colnames(temp)==FALSE|
       "BBP700"%in% colnames(temp)==FALSE){
      print("insufficient data")
      write.csv(temp, file = fileList[i],
                row.names=FALSE)
      next
    }
    
    for(j in 1:length(unique(temp$date))){
      print(j)
      tempFlag = 1
      currentDate_char <- unique(temp$date)[j]
      range = which(temp$date==unique(temp$date)[j])
      currentDate <- as.POSIXct(unique(temp$date)[j])
      currentMLD <- temp$MLD_003[range][1]
      currentProvince <- temp$province[range][1]
      
      if(all(is.na(temp$CHLA[range]))|
         all(is.na(temp$province[range]))|
         is.na(currentMLD)){
        next
      }
      
      temp$MLD_INT_PAR[range] <- avgMldPAR(lon = temp$longitude[range][1],
                                           lat = temp$latitude[range][1],
                                           time_posix = currentDate,
                                           par_profile = temp$DOWNWELLING_PAR[range],
                                           depth_profile = temp$PRES[range],
                                           MLD_depth = currentMLD)
      print(temp$MLD_INT_PAR[range][1])
    }
    
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
    rm(temp)
  }
}

intParFinder <- function(fileList, replace = FALSE){
  for(i in 1:length(fileList)){
    print(paste(i,fileList[i]))
    temp<-read.csv(fileList[i])
    currentFloat = temp$float[1]
    temp$IPAR20 = NA
    
    if(length(grep(currentFloat,parFiles))==0|
       #length(unique(temp$date))<20|
       "CHLA"%in% colnames(temp)==FALSE|
       "BBP700"%in% colnames(temp)==FALSE){
      print("insufficient data")
      write.csv(temp, file = fileList[i],
                row.names=FALSE)
      next
    }
    
    parData<-read.csv(parFiles[grep(currentFloat,parFiles)]) %>% 
      arrange(floatID,date,PRES)
    
    for(j in 1:length(unique(temp$date))){
      print(j)
      
      range = which(temp$date==unique(temp$date)[j])
      parRange = which(as.character(parData$date)==as.character(unique(temp$date)[j]))
      #Set parTemp as estimated par values within parRange (values already interpolated to 1m)
      
      #Set presTemp as estimated PRES values within parRange (values already interpolated to 1m)
      presTemp = parData$PRES[parRange]
      parTemp = parData$DOWNWELLING_PAR[parRange]
      
      currentDate <- as.POSIXct(unique(temp$date)[j])
      currentBasin <- temp$basin[range][1]
      currentMLD <- temp$MLD_003[range][1]
      currentProvince <- temp$province[range][1]
      
      if(all(is.na(parTemp))|
         (is.na(currentProvince))){
        next
      }
      
      tempPARprof <- approx(presTemp,parTemp,1:max(presTemp))$y
      tempPRESprof <-  approx(presTemp,parTemp,1:max(presTemp))$x
      
      #plot(tempPRESprof,tempPARprof)
      
      
      parValues = c()
      output = 1000
      index = 0
      
      while(output > 10){
        index = index + 1
        output = (daily_PAR_given_depth(lon = temp$longitude[range][1],
                                             lat = temp$latitude[range][1],
                                             time_posix = currentDate,
                                             par_profile = tempPARprof,
                                             depth_profile = tempPRESprof,
                                             given_depth = tempPRESprof[index]))
        if(is.na(output)){
          temp$IPAR20[range] = NA
          break
        }
        parValues <- c(parValues,output)
      }
      
      if(all(parValues<20)){
        temp$IPAR20[range] = NA
        next
      }
      
      ipar20 = round(parValues[which(abs(parValues-20)==min(abs(parValues-20),na.rm=TRUE))],2)
      ipar20depth = median(tempPRESprof[which(abs(parValues-20)==min(abs(parValues-20),na.rm=TRUE))])
      print(paste("PAR:",ipar20,"depth:",ipar20depth))
    
     temp$IPAR20[range]=ipar20depth
    }
    write.csv(temp, file = fileList[i],
              row.names=FALSE)
  }
}

parProfile(plot=FALSE)
parFiles<-list.files("data/parData",pattern=".csv",
                     full.names = T)
#Once nutricline parameters and MLD have been calculated, the following functions
#estimate daily integrated PAR at different nutricline depths and within the mixed layer
intPAR(fileList,parFiles,plot=TRUE)
intMldPAR(fileList,replace=TRUE)
#intParFinder(fileList)

