#Script loads profile index, identifies floats w/ CHL, identifies whether
#sufficient data exists for download, and determines whether more recent version
#is available than on disk

basinFinder <- function(lon,lat){
  require(sp)
  require(readxl)
  #browser()
  if(length(grep("regionPolygonsInclusive",list.files("supportingFiles/")))==0){
    print("Basin defitions file not found in supporting files directory")
    return()
  }
  basinGen <- read_excel("supportingFiles/regionPolygonsInclusive.xlsx",
                         sheet="gen")
  bgNames<- unique(basinGen$area)
  basins = rep("BadCoords",length(lat))
  for(i in 1:length(lat)){
    if(floor(i/1e4) == i/1e4){
      print(paste(round(i/length(lat),2)*100,"% complete"))
    }
    for(j in 1:length(bgNames)){
      bgCoords <- subset(basinGen, area == bgNames[j])
      if(point.in.polygon(lon[i],lat[i],bgCoords$lon,bgCoords$lat) >= 1){
        if(bgNames[j] %in% c("SEPac","SWPac")){
          basins[i] = "Pacific_S"
          break
        }
        if(bgNames[j] %in% c("NEPac","NWPac")){
          basins[i] = "Pacific_N"
          next
        }
        if(bgNames[j] %in% c("EqWPac","EqEPac")){
          basins[i] = "Pacific_Eq"
          break
        }
        else{
          basins[i] = bgNames[j]
          break
        }
      }
    }
  }
  return(basins)
}

###########################################################################

bathGet <- function(lat,lon){
  if(length(grep("GEBCO",list.files("supportingFiles/")))==0){
    cat(paste("GEBCO bathymetry file not found in supporting files directory.",
                "Download from https://www.gebco.net/data_and_products/gridded_bathymetry_data/",sep="\n"))
    return()
  }

  bath<- read.nc(open.nc("supportingFiles/GEBCO_2020.nc"))
  if(length(lat) != length(lon)){
    print("Coordinate vectors of differing lengths")
    return()
  }
  elev = rep(NA,length(lat))
  for(i in 1:length(lat)){
    if(floor(i/1e4) == i/1e4){
      print(paste(round(i/length(lat),2)*100,"% complete"))
    }
    if(is.na(lat)[i] | is.na(lon)[i]){
      next}
    bathLatIndex <- abs(bath$lat-lat[i])
    bathLonIndex <- abs(bath$lon-lon[i])
    bathLatIndex <- min(which(bathLatIndex == min(bathLatIndex)))
    bathLonIndex <- min(which(bathLonIndex == min(bathLonIndex)))
    #longitude is stored in rows, latitude stored in columns
    elev[i] <- bath$elevation[bathLonIndex,bathLatIndex]
  }
  return(elev)
}

############################################################################

vectorize<-function(input,varName, nP,nD){
  #browser()
  require(stringr)
  if(grepl("^PROFILE",varName)==TRUE){
    outputVector = rep(unlist(str_split(input,""),nD),each=nD)
    return(as.character(outputVector))
  }
  
  if(grepl("QC$",varName)==TRUE & 
     grepl("^PROFILE",varName)!=TRUE){
    outputVector<-rep(NA,nP*nD)
    for(i in 1:length(input)){
      if(is.na(nchar(input[i],allowNA = TRUE))){
        outputVector[(((i-1)*nD)+1):(i*nD)] = NA
      }
      else if(nchar(input[i],allowNA = TRUE) == nD){
        outputVector[(((i-1)*nD)+1):(i*nD)] = unlist(str_split(input[i],""))
      }
      else{
        outputVector[(((i-1)*nD)+1):(i*nD)] = NA
      }
    }
    outputVector[which(outputVector==" ")] <- NA
    return(as.numeric(outputVector))
  }
  
  if(length(input)==(nP*nD)){
    return(as.numeric(as.vector(input)))
  }
}

############################################################################

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


##############################################################################


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

###############################################################################

isolumes <- function(fileList,parFiles,
                     inputDirectory = "data/finalCSV/",
                     outputDirectory = "data/finalCSV/"){
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read_csv(paste(inputDirectory,fileList[i],sep=""))
    currentFloat = temp$float[1]
    currentBasin = temp$basin[1]
    
    temp$PAR_FLAG = NA
    temp$ISO_10 = NA
    temp$ISO_1 = NA
    temp$ISO_01 = NA
    
    if("ISO_01_ORIG" %in% colnames(temp)){
      temp<-temp %>% 
        select(-ISO_01_ORIG)
    }
    
    #If the current float has a corresponding PAR file, open PAR file as parData
    if(length(grep(currentFloat,parFiles))==1){
      parData<-read_csv(parFiles[grep(currentFloat,parFiles)]) %>% 
        arrange(floatID,date,PRES)
    }
    
    #If the current float does not have a corresponding PAR file, skip. Save so
    #that all files have the isolume columns, though.
    if(length(grep(currentFloat,parFiles))!=1){
      write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
                row.names = FALSE)
      next
    }
    for(j in 1:length(unique(temp$date))){
      print(j)
      #Identify current date
      currentDate = unique(temp$date)[j]
      #Identify rows in data in temp matching current date
      range = which(temp$date==currentDate)
      #Identify rows in parData matching current date
      parRange = which(parData$date==currentDate)
      #Set parTemp as estimated par values within parRange (values already interpolated to 1m)
      parTemp = parData$DOWNWELLING_PAR[parRange]
      
      parFlag = parData$PAR_FLAG[parRange][1]
      #Set presTemp as estimated PRES values within parRange (values already interpolated to 1m)
      presTemp = parData$PRES[parRange]
      
      #If all values in parTemp are missing, or if all values are 0, skip profile
      if(all(is.na(parTemp))|
         all(parTemp==0,na.rm=TRUE)){
        next
      }
      
      #Calculate depths where PAr is 10 percent, 1 percent, 0.1 percent of PAR at 1 meter
      iso10profile = (parTemp/parTemp[1]) - 0.1
      iso1profile = (parTemp/parTemp[1]) - 0.01
      iso01profile = (parTemp/parTemp[1]) - 0.001
      
      temp$ISO_10[range] = presTemp[iso10profile<=0][1]
      temp$ISO_1[range] = presTemp[iso1profile<=0][1]
      temp$ISO_01[range] = presTemp[iso01profile<=0][1]
      temp$PAR_FLAG[range] = parFlag
    }
    write_csv(temp, paste(outputDirectory,fileList[i],sep=""))
  }
}


################################################################################

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

################################################################################

oxySat <- function(doxy,temp, sal){
  #Doxy: umol kg-1
  #Calculates oxygen solubility (umol kg-1, based on equations in Garcia + Gordon 1992),
  #apparent oxygen utilization, and saturation for given values of doxy, temp and salinity
  A0 = 5.80818
  A1 = 3.20684
  A2 = 4.11890
  A3 = 4.93845
  A4 = 1.01567
  A5 = 1.41575
  B0 = -7.01211e-3
  B1 = -7.25958e-3
  B2 = -7.93334e-3
  B3 = -5.54491e-3
  C0 = -1.32412e-7
  
  temp.adj <- log((298.15-temp)/(273.15+temp))
  logSol <- A0 + (A1*temp.adj) + (A2 * (temp.adj^2)) +
    (A3 * (temp.adj^3)) + (A4 * (temp.adj^4)) + (A5 * (temp.adj^5)) +
    (sal * (B0 + (B1*temp.adj) + (B2*(temp.adj^2)) + (B3 * (temp.adj^3)))) +
    (C0 * (sal^2))
  sol = exp(logSol)
  aou = sol - doxy
  sat = (doxy/sol)*100
  
  #output<-data.frame(sol,aou,sat)
  #return(output)
  return(aou)
}



