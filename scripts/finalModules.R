require(gsw)
require(tidyverse)
require(lubridate)
require(zoo)
require(cowplot)

###############################################################################
###############################################################################
###################SCRIPTS SUPPORTING PARAMETER DETERMINATION##################
###############################################################################
###############################################################################
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

bathGet <- function(lat,lon,bath){
  bathLatIndex <- abs(bath$lat-lat)
  bathLonIndex <- abs(bath$lon-lon)
  bathLatIndex <- min(which(bathLatIndex == min(bathLatIndex)))
  bathLonIndex <- min(which(bathLonIndex == min(bathLonIndex)))
  #longitude is stored in rows, latitude stored in columns
  depth <- bath$elevation[bathLonIndex,bathLatIndex]
  return(depth)
}

ncliner<-function(NUTE,PRES,THRESH,FACTOR){
  #browser()
  NUTE <- rollmedian(NUTE,7,"extend")
  presTemp <- approx(PRES,NUTE,1:max(PRES,na.rm=TRUE))$x
  nuteTemp <-  na.fill(approx(PRES,NUTE,1:max(PRES,na.rm=TRUE))$y,"extend")

  #Calculate nitracline parameters based on CanyonB estimates
  #Calculate surface nutrient as the mean of nutrient measurements within first 10 meters
  nuteSurf = mean(nuteTemp[presTemp<=5],na.rm=TRUE)
  #Calculate difference between nutrient at each depth and surface measurements; 
  #Nutricline defined as depth where nutrient is 1 uM greater than at the surface
  #Nutricline (calculated following Cornec 2020) is derived from region between 1 and 1.25 times the nutricline depth; 
  #here, region between 1 and 1.25 times the nutricline index is used, to ensure sufficient data for regressions at shallow depths 
  #Calculate difference between nutrient concentration, surface nutrients, and threshold (1 uM).
  #Value closest to 0 corresponds to the nutricline
  nuteDiff = (nuteTemp-nuteSurf)
  nutriclineDepth<-presTemp[nuteDiff >= THRESH][1]
  
  #Determine vector index corresponding to nitracline depth
  nutriclineIndex_1 = which(presTemp==nutriclineDepth)
  #In some cases nutrients may not exceed the threshold, in this case, set nutricline depth
  #to NA (because there is no nutricline), and slope to 0
  if(is.na(nutriclineDepth)){
    nutriclineDepth = NA
    nutriclineSlope = 0
    nutriclineInt = 0
    nutriclineR2=NA
    return(c(nutriclineDepth,nutriclineSlope,nutriclineInt,nutriclineR2))
    }
  #Nutricline depth is the depth corresponding to the minimum value of nuteDiff
  #Note that max() is used to ensure only one value is returned in cases where there are duplicate
  #values for nuteDiff.
  #Determine vector index corresponding to the prduct of nitraclineIndex_1 and the desired factor (1.25)
  # nutriclineIndex_2 = nutriclineIndex_1 + 10
  #nutriclineIndex_2 = ceiling(nutriclineIndex_1*FACTOR)
  nutriclineIndex_2 = min(which(abs(presTemp-(nutriclineDepth*1.25))==min(abs(presTemp-(nutriclineDepth*1.25)))))
  # 
  if(nutriclineIndex_2 - nutriclineIndex_1 < 3){
    nutriclineIndex_2 = nutriclineIndex_1 + 3
  }
  #Perform rolling regression 3 depths above and below 
  regRange <- seq(-5,5,1)
  #In homongenous water columns, nitraclineIndex_1 may be extremely deep; in these instances,
  #set nitraclineIndex_2 to the last index value in the profile
  #If nitraclineIndex is too large, it won't be possible to perform rolling regression loop
  #below (because nitraclineIndex_2 is set to the length of the depth vector); in this case,
  #set range to 0 so that the regression loop only performs one regression between 
  #nitraclineIndex_1 and nitraclineIndex2
  if(nutriclineIndex_1 > length(presTemp)*.8){
    nutriclineIndex_2 = length(presTemp)
    regRange = 0
  }
  
  #If nutriclineIndex_1 is too deep to allow for regression to be performed, set slope to 0 
  #and other values to NA
  if(nutriclineIndex_1 > length(presTemp)-3){
    nutriclineDepth = presTemp[nutriclineIndex_1]
    nutriclineSlope = 0
    nutriclineInt = NA
    nutriclineR2 = NA
    return(c(nutriclineDepth,nutriclineSlope,nutriclineInt,nutriclineR2))
  }
  
  if(nutriclineIndex_1 < 3){
    regRange <- seq(-nutriclineIndex_1+1,nutriclineIndex_1-1,1)
  }
  
  #Subset PRES and NUTE based on nutriclineIndex1 & nutriclineIndex2
  #presSubset <- PRES[nutriclineIndex_1:nutriclineIndex_2]
  #nuteSubset <-NUTE[nutriclineIndex_1:nutriclineIndex_2]
  
  slopeList <- c()
  rsqList <- c()
  intList<-c()
  for(i in regRange){
    presSubset <- presTemp[(nutriclineIndex_1+i):(nutriclineIndex_2+i)]
    nuteSubset <-nuteTemp[(nutriclineIndex_1+i):(nutriclineIndex_2+i)]
    #print(nuteSubset)
    temp<-summary(lm(nuteSubset~presSubset))
    slope <- temp$coefficients[2,1]
    int<-temp$coefficients[1,1]
    rsq <- temp$r.squared
    slopeList <- c(slopeList,slope)
    rsqList <- c(rsqList,rsq)
    intList<-c(intList,int)
  }
  
  #Perform regressions comparing PRES and NUTE; calculate mean in case there
  #are multiple identical slope values within slope list
  nutriclineSlope = mean(slopeList[which(slopeList == max(slopeList))],na.rm=TRUE)
  nutriclineInt = mean(intList[which(slopeList == max(slopeList))],na.rm=TRUE)
  nutriclineR2 = mean(rsqList[which(slopeList == max(slopeList))],na.rm=TRUE)
  
  return(c(nutriclineDepth,nutriclineSlope,nutriclineInt,nutriclineR2))
}

match_up_sat<-function(month,year,lat,lon, path_sat_files) {
  #browser()
  # Inputs : 
  #   month : month of the profile (value in format MM)
  #   year : year of the profile ( value in format YYYY)
  #   lon : longitude of the profile
  #   lat : latitude of the profile
  #   path_sat_file: path to folder where the .nc files are downloaded
  #   
  #   Output :
  #     chla_sat_value: closest geographical chla sat value with the profile
  
  library(ncdf4)
  library(stringr)
  #add a 0 if the month has less than 2 charactters
  month<-str_pad(month,2, pad = "0")
  
  #empty sat variable
  sat<-NA
  
  filename = paste(path_sat_files,"ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-",
                   year,month,"-fv5.0.nc",sep="")
  # 
  # if((filename %in% list.files(path_sat_files))==FALSE){
  #   return(NA)
  # }
  # 
  # open sat file
  sat<-nc_open(filename,readunlim=FALSE,write=FALSE)
  
  
  # read chla, lon , and lat from sat
  chl_sat<-NA
  lat_sat<-NA
  lon_sat<-NA
  
  chl_sat<-ncvar_get(sat,"chlor_a")
  lat_sat<-ncvar_get(sat,"lat")
  lon_sat<-ncvar_get(sat,"lon")
  
  #close the nc file
  nc_close(sat)
  
  coord_lat<-NA
  coord_lon<-NA
  
  # find the closest sat point from the geolocation of the profile
  coord_lat<-which(abs(lat_sat-lat)==
                     min(abs(lat_sat-lat),na.rm=T))[1]
  coord_lon<-which(abs(lon_sat-lon)==
                     min(abs(lon_sat-lon),na.rm=T))[1]
  
  chla_sat_value<-NA
  chla_sat_mean<-NA
  chla_sat_sd<-NA
  
  chla_sat_value<-chl_sat[coord_lon,coord_lat]
  
  coord_lat = c((coord_lat-1):(coord_lat+1))
  coord_lon = c((coord_lon-1):(coord_lon+1))
  
  chla_sat_mean<-mean(chl_sat[coord_lon,coord_lat],na.rm=TRUE)
  chla_sat_sd<-sd(chl_sat[coord_lon,coord_lat],na.rm=TRUE)
  
  return(c(chla_sat_value,
           chla_sat_mean,
           chla_sat_sd))
}

integrator<-function(param,depths,start,stop){
  require(zoo)
  #browser()
  
  if(sum(is.na(depths)==FALSE)<=2|
     sum(is.na(param)==FALSE)<=2|
     is.na(start)|
     is.na(stop)){
    return(NA)
  }
  
  depthTemp <- approx(depths,param,start:stop)$x
  paramTemp <-  approx(depths,param,start:stop)$y
  
  param<-na.fill(paramTemp,"extend")
  

  #Selects parameter measurements greater than and including minDepth, less than maxDepth;
  #this ignores the maximum depth in any given profile, but this is presumably less consequential
  #than ignoring the first value (or including a given value twice if both ends are inclusive)
  #param = paramTemp[depthTemp>=start & depthTemp<=stop]
  
  return(sum(param,na.rm=TRUE))
}

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
###############################################################################
########################FUNCTIONS TO DETERMINE PARAMETERS######################
###############################################################################
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


##################################MLD#########################################
mld<-function(fileList,plot=FALSE,
              inputDirectory = "data/Chla_BBP_data/finalCSV/",
              outputDirectory = "data/Chla_BBP_data/finalCSV/"){
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
############################MLD INTEGRATED NUTRIENTS############################
intNewts<-function(fileList,
                   inputDirectory = "data/finalCSV/",
                   outputDirectory = "data/finalCSV/"){
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep=""))%>% 
      arrange(date,PRES)
    
    temp$AVG_MLD_NO3 = NA
    temp$AVG_MLD_CAN_NO3 = NA
    temp$AVG_MLD_CAN_PO4 = NA
    temp$AVG_MLD_CAN_SiO4 = NA
    
    if("AVG_MLD_N03" %in% colnames(temp)){
      temp<-temp %>% 
        select(-AVG_MLD_N03)
    }
    
    for(j in 1:length(unique(temp$date))){
      print(j)

      range = temp$date==unique(temp$date)[j]
      currentMLD = mean(temp$MLD_003[range],1)
      if(is.na(currentMLD)|
         currentMLD==1){
        print("MLD too shallow")
        next
      }
      currentSig <- ((temp$SIG[range]+1000)/1000)
      temp$NITRATE[range][temp$NITRATE[range]<0]=0
      temp$CAN_NO3[range][temp$CAN_NO3[range]<0]=0
      temp$CAN_PO4[range][temp$CAN_PO4[range]<0]=0
      temp$CAN_SiO4[range][temp$CAN_SiO4[range]<0]=0
      
      temp$AVG_MLD_NO3[range] =integrator(param = temp$NITRATE[range]*currentSig,depths = temp$PRES[range],
                                      start = 1,stop =currentMLD)/currentMLD
      temp$AVG_MLD_CAN_NO3[range] = integrator(param = temp$CAN_NO3[range]*currentSig,depths=temp$PRES[range],
                                           start = 1,stop =currentMLD)/currentMLD
      temp$AVG_MLD_CAN_PO4[range] = integrator(param = temp$CAN_PO4[range]*currentSig,depths=temp$PRES[range],
                                           start = 1,stop =currentMLD)/currentMLD
      temp$AVG_MLD_CAN_SiO4[range] = integrator(param = temp$CAN_SiO4[range]*currentSig,depths=temp$PRES[range],
                                            start = 1,stop =currentMLD)/currentMLD

        
    }
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }}

bbpInt<-function(fileList,
                 inputDirectory = "data/finalCSV/",
                 outputDirectory = "data/finalCSV/"){
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep=""))
    
    currentFloat = temp$float[1]
    temp$INT_BBP_MESO_0_1 = NA
    temp$INT_BBP_MESO_1_2 = NA
    temp$INT_BBP_MESO_2_3 = NA
    temp$INT_BBP_MESO_3_4 = NA
    temp$INT_BBP_MESO_0_4 = NA
    temp$AVG_BBP_MESO_0_1 = NA
    temp$AVG_BBP_MESO_1_2 = NA
    temp$AVG_BBP_MESO_2_3 = NA
    temp$AVG_BBP_MESO_3_4 = NA
    temp$AVG_BBP_MESO_0_4 = NA
    temp$AVG_BBP_MESO_SLOPE = NA
    temp$AVG_BBP_MESO_SLOPE_R = NA
    temp$INT_BBP_PZ = NA
    temp$AVG_BBP_PZ = NA
    temp$INT_CHL_PZ = NA
    temp$AVG_CHL_PZ = NA

    if("INT_BBP_PZ_YR" %in% colnames(temp)){
      temp<-temp %>% 
        select(-INT_BBP_PZ_YR,
               -INT_BBP_MESO_YR)
    }
    
    text <- "PAR data,"
    if(length(grep(currentFloat,parFiles))==0){
      text <- "No PAR data"
    }
    
    for(j in 1:length(unique(temp$date))){
      print(j)
      currentDate = unique(temp$date)[j]
      range = which(temp$date==currentDate)
      if(all(is.na(temp$ISO_01[range]))|
         all(is.na(temp$BBP700[range]))|
         all(temp$BBP700[range]==0)|
         all(temp$ISO_01[range]==0)){
        print(paste(text,"no ISO_01"))
        next
      }
      currentFloat = temp$floatID[range][1]
      currentMLD = temp$MLD_003[range][1]
      currentISO_01 = temp$ISO_01[range][1]
      currentISO_1 = temp$ISO_1[range][1]
      currentISO_10 = temp$ISO_10[range][1]
      currentPz = max(currentMLD,currentISO_01,na.rm=TRUE)
      currentMax = max(temp$PRES[range],na.rm=TRUE)
      mesoDepths = seq(from = currentPz,
                       to = currentMax,
                       by = (currentMax-currentPz)/(4))

      temp$INT_BBP_MESO_0_1[range] = integrator(temp$BBP700[range],
                                                temp$PRES[range],
                                                start = currentPz,
                                                stop=mesoDepths[2])
    
      temp$AVG_BBP_MESO_0_1[range] = temp$INT_BBP_MESO_0_1[range]/(mesoDepths[2]-currentPz)
      
      temp$INT_BBP_MESO_1_2[range] = integrator(temp$BBP700[range],
                                                temp$PRES[range],
                                                start = mesoDepths[2],
                                                stop=mesoDepths[3])
      
      temp$AVG_BBP_MESO_1_2[range] = temp$INT_BBP_MESO_1_2[range]/(mesoDepths[3]-mesoDepths[2])
      
      temp$INT_BBP_MESO_2_3[range] = integrator(temp$BBP700[range],
                                                temp$PRES[range],
                                                start =mesoDepths[3],
                                                stop=mesoDepths[4])
      
      temp$AVG_BBP_MESO_2_3[range] = temp$INT_BBP_MESO_2_3[range]/(mesoDepths[4]-mesoDepths[3])
      
      temp$INT_BBP_MESO_3_4[range] = integrator(temp$BBP700[range],
                                                temp$PRES[range],
                                                start = mesoDepths[4],
                                                stop=mesoDepths[5])
      
      temp$AVG_BBP_MESO_3_4[range] = temp$INT_BBP_MESO_3_4[range]/(mesoDepths[5]-mesoDepths[4])
      
      temp$INT_BBP_MESO_0_4[range] = integrator(temp$BBP700[range],
                                                temp$PRES[range],
                                                start = currentPz,
                                                stop=mesoDepths[5])
      
      temp$AVG_BBP_MESO_0_4[range] = temp$INT_BBP_MESO_0_4[range]/(mesoDepths[5]-currentPz)
      
      temp$INT_BBP_PZ[range] = integrator(temp$BBP700[range],
                                          temp$PRES[range],
                                          start = 1,
                                          stop=currentPz)
      
      temp$AVG_BBP_PZ[range] = temp$INT_BBP_PZ[range]/currentPz
      
      temp$INT_CHL_PZ[range] = integrator(temp$CHLA[range],
                                          temp$PRES[range],
                                          start = 1,
                                          stop=currentPz)
      
      temp$AVG_CHL_PZ[range] = temp$INT_CHL_PZ[range]/currentPz

      bbpValues=c(temp$AVG_BBP_MESO_0_1[range][1],temp$AVG_BBP_MESO_1_2[range][1],
                  temp$AVG_BBP_MESO_2_3[range][1],temp$AVG_BBP_MESO_3_4[range][1]) * 1000
      
      if(length(unique(bbpValues)) < 3){
        next}
      
      depths = rollmean(mesoDepths, k = 2)
      
      temp$AVG_BBP_MESO_SLOPE[range]= summary(lm(bbpValues~depths))$coefficients[2,1]
      temp$AVG_BBP_MESO_SLOPE_R[range]= summary(lm(bbpValues~depths))$adj.r.squared
    }
    
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }
}

canyonB <- function(fileList,replace=FALSE,
                    inputDirectory = "data/finalCSV/",
                    outputDirectory = "data/finalCSV/"){
  #browser()
  Sys.setenv(tz = "UTC")
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep="")) %>% 
      filter(is.na(date)==FALSE)
    
    if("CHLA"%in%colnames(temp)==FALSE | "BBP700" %in% colnames(temp)==FALSE){
      temp$CAN_NO3 = NA
      temp$CAN_NO3_CI = NA
      temp$CAN_PO4 = NA
      temp$CAN_PO4_CI = NA
      temp$CAN_SiO4 = NA
      temp$CAN_SiO4_CI = NA
      write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
                row.names = FALSE)
      next
    }

    temp$CAN_NO3 = NA
    temp$CAN_NO3_CI = NA
    temp$CAN_PO4 = NA
    temp$CAN_PO4_CI = NA
    temp$CAN_SiO4 = NA
    temp$CAN_SiO4_CI = NA
    for(j in 1:length(unique(temp$date))){
      print(j)
      currentDate = unique(temp$date)[j]
      range = temp$date==currentDate
      if(all(is.na(temp$DOXY[range])) | all(is.na(temp$PRES[range])) |
         all(is.na(temp$PSAL[range])) | all(is.na(temp$TEMP[range]))){
        next}
      
      currentBasin<-temp$basin[range][1]
      currentLat <- temp$latitude[range][1]
      currentLon <- temp$longitude[range][1]
      if(currentBasin %in% c("Med_E","Med_W")==TRUE){
        temp$CAN_NO3[range] = CANYON_MED_NO3_v4(date = currentDate,
                                                lat = currentLat,
                                                lon = currentLon,
                                                pres = temp$PRES[range],
                                                temp = temp$TEMP[range],
                                                psal = temp$PSAL[range],
                                                doxy = temp$DOXY[range])
        
        temp$CAN_PO4[range] = CANYON_MED_PO4_v4(date = currentDate,
                                                lat = currentLat,
                                                lon = currentLon,
                                                pres = temp$PRES[range],
                                                temp = temp$TEMP[range],
                                                psal = temp$PSAL[range],
                                                doxy = temp$DOXY[range])
        
        temp$CAN_SiO4[range] = CANYON_MED_SiOH4_v4(date = currentDate,
                                                   lat = currentLat,
                                                   lon = currentLon,
                                                   pres = temp$PRES[range],
                                                   temp = temp$TEMP[range],
                                                   psal = temp$PSAL[range],
                                                   doxy = temp$DOXY[range])
      }
      if(currentBasin %in% c("Med_E","Med_W")==FALSE){
        canyonB = CANYONB(date = currentDate,
                          lat = currentLat,
                          lon = currentLon,
                          pres = temp$PRES[range],
                          temp = temp$TEMP[range],
                          psal = temp$PSAL[range],
                          doxy = temp$DOXY[range],
                          param=c('NO3','PO4','SiOH4'))
        temp$CAN_NO3[range]=canyonB$NO3
        temp$CAN_NO3_CI[range]=canyonB$NO3_ci
        temp$CAN_PO4[range]=canyonB$PO4
        temp$CAN_PO4_CI[range]=canyonB$PO4_ci
        temp$CAN_SiO4[range] = canyonB$SiOH4
        temp$CAN_SiO4_CI[range] = canyonB$SiOH4_ci
      }
    }
    write.csv(temp,paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }
}

chlSat<-function(fileList,replace=FALSE,
                 inputDirectory = "data/finalCSV_params/",
                 outputDirectory = "data/finalCSV_params/"){
  
  satPath<-"/Users/nicholasbock/Projects/argo/data/satChl/"
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep=""))
    
    temp$CHL_SAT = NA
    for(j in 1:length(unique(temp$date))){
      print(j)
      currentDate = unique(temp$date)[j]
      if(is.na(currentDate)|year(currentDate)==2021){
        next
      }
      currentLat = mean(temp$latitude[temp$date==currentDate],na.rm=TRUE)
      currentLon = mean(temp$longitude[temp$date==currentDate],na.rm=TRUE)
      satData <- match_up_sat(month = month(currentDate),
                              year = year(currentDate),
                              lat = currentLat,
                              lon = currentLon,
                              path_sat_files = satPath)

      temp$CHL_SAT[temp$date==currentDate] <- satData[1]
      temp$CHL_SAT_AVG[temp$date==currentDate] <- satData[2]
      temp$CHL_SAT_SD[temp$date==currentDate] <- satData[3]
      }

    write.csv(temp,paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }   
}

nclines<-function(fileList,plot=FALSE,
                  inputDirectory = "data/finalCSV_params/",
                  outputDirectory = "data/finalCSV_params/"){
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep="")) %>% 
      mutate(PRES = as.numeric(PRES)) %>% 
      arrange(jDay,PRES)
    currentBasin<-temp$basin[1]
    currentFloat<-temp$floatID[1]
    
    if("NITRATE" %in% colnames(temp)==FALSE){
      temp$NITRATE = NA
    }
    
    temp$NCLINE_DEPTH = NA
    temp$NCLINE_SLOPE = NA
    temp$NCLINE_R2 = NA
    
    temp$CAN_NCLINE_DEPTH = NA
    temp$CAN_NCLINE_SLOPE = NA
    temp$CAN_NCLINE_R2 = NA
    
    temp$CAN_SCLINE_DEPTH = NA
    temp$CAN_SCLINE_SLOPE = NA
    temp$CAN_SCLINE_R2 = NA
    
    temp$CAN_PCLINE_DEPTH = NA
    temp$CAN_PCLINE_SLOPE = NA
    temp$CAN_PCLINE_R2 = NA
      
    for(j in 1:length(unique(temp$date))){
      print(j)
      currentDate = unique(temp$date)[j]
      currentFloat = temp$floatID[temp$date==currentDate][1]
      range = temp$date==currentDate 
      noFactor = 1

      if(all(is.na(temp$CAN_NO3[range]))==FALSE){
        siFactor <- round(max(temp$CAN_SiO4[range],na.rm=TRUE)/max(temp$CAN_NO3[range],na.rm=TRUE),3)
        pFactor <- round(max(temp$CAN_PO4[range],na.rm=TRUE)/max(temp$CAN_NO3[range],na.rm=TRUE),3)
        #NCliner takes nitrate and pressure as arguments, and returns:
        #NCline depth, slope, int, and R2
        can_nclineParams<-ncliner(temp$CAN_NO3[range],temp$PRES[range],noFactor,1.25)
        sclineParams<-ncliner(temp$CAN_SiO4[range],temp$PRES[range],siFactor,1.25)
        pclineParams<-ncliner(temp$CAN_PO4[range],temp$PRES[range],pFactor,1.25)

        temp$CAN_NCLINE_DEPTH[range] = can_nclineParams[1]
        temp$CAN_NCLINE_SLOPE[range] = can_nclineParams[2]
        temp$CAN_NCLINE_R2[range] = can_nclineParams[4]

        temp$CAN_SCLINE_DEPTH[range] = sclineParams[1]
        temp$CAN_SCLINE_SLOPE[range] = sclineParams[2]
        temp$CAN_SCLINE_R2[range] = sclineParams[4]

        temp$CAN_PCLINE_DEPTH[range] = pclineParams[1]
        temp$CAN_PCLINE_SLOPE[range] = pclineParams[2]
        temp$CAN_PCLINE_R2[range] = pclineParams[4]
      }
      
      if(all(is.na(temp$NITRATE[range]))==FALSE){
        nclineParams<-ncliner(temp$NITRATE[range],temp$PRES[range],1,1.25)
        temp$NCLINE_DEPTH[range] = nclineParams[1]
        temp$NCLINE_SLOPE[range] = nclineParams[2]
        temp$NCLINE_R2[range] = nclineParams[4]
      }
        
      if(plot==TRUE & all(is.na(temp$CAN_NO3[range]))==FALSE){
        PO4plot <- ggplot() +
          geom_line(aes(x = temp$PRES[range], y = temp$CAN_PO4[range])) +
          coord_flip() +
          geom_vline(xintercept = pclineParams[1], linetype =2) +
          geom_abline(intercept=pclineParams[3],slope=-pclineParams[2],col="red")+
          scale_x_reverse() +
          labs(y = "PO4", x = "depth") +
          ggtitle(paste(paste("R =",round(pclineParams[4],2)," || Factor = ",pFactor,sep=""),
                        paste("Slope =", round(pclineParams[2],2)),sep="\n")) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 8))
        
        NO3plot <- ggplot() +
          geom_line(aes(x = temp$PRES[range], y = temp$CAN_NO3[range])) +
          coord_flip() +
          geom_vline(xintercept = can_nclineParams[1], linetype =2) +
          geom_abline(intercept=can_nclineParams[3],slope=-can_nclineParams[2],col="red")+
          scale_x_reverse() +
          ggtitle(paste(paste("R =",round(can_nclineParams[4],2)," || Factor = ",noFactor,sep=""),
                        paste("Slope =", round(can_nclineParams[2],2)),sep="\n")) +
          labs(y = "NO3", x = "") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 8))
        
        SiOplot <- ggplot() +
          geom_line(aes(x = temp$PRES[range], y = temp$CAN_SiO4[range])) +
          coord_flip() +
          geom_vline(xintercept = sclineParams[1], linetype =2) +
          geom_abline(intercept=sclineParams[3],slope=-sclineParams[2],col="red")+
          scale_x_reverse() +
          labs(y = "SiO4", x = "") +
          ggtitle(paste(paste("R =",round(sclineParams[4],2)," || Factor = ",siFactor,sep=""),
                        paste("Slope =", round(sclineParams[2],2)),sep="\n")) +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 8))
        
        outputPlot <- plot_grid(PO4plot,NO3plot,SiOplot,
                                ncol = 3)
        
        ggsave(paste("figures/parameterTestPlots/nclineSlopeTest/",
                     currentCluster,"_",currentFloat,"_",
                     str_pad(j,3,side="left",pad="0"),".png",sep=""),
               plot=outputPlot,dpi = 150)
    }
  }
    
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }
}

dcmDepth<-function(fileList,plot=FALSE,
                   inputDirectory = "data/finalCSV_params/",
                   outputDirectory = "data/finalCSV_params/"){
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep="")) %>% 
      arrange(cycleNumber,PRES)
    
    if("BBP_MED_MAX" %in% colnames(temp)){
      temp<-temp[,-which(colnames(temp)=="BBP_MED_MAX")]
      temp<-temp[,-which(colnames(temp)=="BBP_MESO_MIN_MAX")]
    }
    
    if("BBP_MAX_MED" %in% colnames(temp)){
      temp<-temp[,-which(colnames(temp)=="BBP_MAX_MED")]
    }
    
    temp$DCM = NA
    temp$DCM_DEPTH = NA
    temp$DCM_MED_MAX = NA
    temp$DCM_DIFF = NA
    temp$DCM_BBP = NA
    temp$DCM_BBP_RATIO = NA
    temp$BBP_MAX = NA
    temp$BBP_MAX_DEPTH =NA
    temp$BBP_MAX_MIN = NA
    temp$BBP_MESO_MAX = NA
    temp$BBP_MESO_MAX_DEPTH = NA
    temp$BBP_MESO_MAX_MIN = NA
    
    if("CHLA_ORIG" %in% colnames(temp)==FALSE|
      all(is.na(temp$CHLA_ORIG))){
      write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
                row.names = FALSE)
      print("insufficient data")
      next
    }
    
    for(j in 1:length(unique(temp$date))){
      print(j)
      isDCM = FALSE
      thresh = 0.05
      
      currentDate = unique(temp$date)[j]
      range = which(temp$date==currentDate)
      currentLat<-temp$latitude[range][1]
      currentFloatID <- temp$floatID[range][1]
      currentBasin<-temp$basin[range][1]
      currentFloat<-temp$floatID[range][1]
      
      plotMaxChl = max(temp$CHLA[temp$floatID==currentFloatID],na.rm=TRUE) * 1.1
      plotMaxBBP = max(temp$BBP700[temp$floatID==currentFloatID],na.rm=TRUE) * 1.1
      
      if(all(is.na(temp$CHLA[range])) |
         min(temp$PRES[range][is.na(temp$CHLA_ORIG[range])==FALSE],na.rm=TRUE)>15){
        print("Missing surface Chl measurements, maybe entire profile")
        next}

      chlProfile<-temp$CHLA_ORIG[range]
      chlProfile[chlProfile<0]=0
      presProfile<-temp$PRES[range][is.na(chlProfile)==FALSE]
      chlProfile<-chlProfile[is.na(chlProfile)==FALSE]
      chlSmoothProfile<-rollmedian(chlProfile,3,align="center","extend")
      
      chlPresTemp<-approx(presProfile,chlProfile,
                       1:max(presProfile))$x
      chlTemp<-approx(presProfile,chlProfile,
                      1:max(presProfile))$y
      chlSmoothTemp<-approx(presProfile,chlSmoothProfile,
                            1:max(presProfile))$y
      
      cmax_index <- which.max(chlSmoothTemp)
      cmax = chlTemp[cmax_index]
      cmed = median(chlTemp[1:10],na.rm=TRUE)
      cdiff = cmax-cmed
      
      temp$DCM[range] = cmax
      temp$DCM_DIFF[range] = cdiff
      temp$DCM_DEPTH[range] = cmax_index
      temp$DCM_MED_MAX[range] = round(cmax/cmed,2)
      
      chlProfilePlot=ggplot()+
        scale_y_continuous(limits=c(0,plotMaxChl))+
        coord_flip() + 
        scale_x_reverse(limits=c(750,0)) + 
        geom_line(aes(x=chlPresTemp, y=chlSmoothTemp),linetype=2,alpha = 0.5) + 
        geom_line(aes(x=chlPresTemp, y=chlTemp)) + 
        theme_bw() + 
        labs(x = "Depth (m)", y = bquote("Chl ("*mg~m^-3*")")) +
        theme(plot.title = element_text(hjust=0.5))
      
      if(round(cmax/cmed,2) >= 1.3 & 
         cdiff > thresh &
         cmax_index < 250){
        isDCM = TRUE
        chlProfilePlot <- chlProfilePlot + 
          geom_vline(xintercept = cmax_index,col ="#91A737")
        }
        
      bbpProfilePlot = NA

      if("BBP700" %in% colnames(temp) & all(is.na(temp$BBP700[range]))==FALSE &
         all(temp$BBP700[range]==0)==FALSE &
         #This next line checks if any measurements are available above 15m
         min(temp$PRES[range][is.na(temp$BBP700_ORIG[range])==FALSE],na.rm=TRUE)<15){
        
        bppProfile<-temp$BBP700_ORIG[range]
        #It is necessary to recreate presProfile because different CHLA and BBP
        #measurements may be missing, which would result in profiles of different lengths
        presProfile<-temp$PRES[range][is.na(bppProfile)==FALSE]
        bbpProfile<-bppProfile[is.na(bppProfile)==FALSE]
        bbpSmoothProfile<-rollmean(rollmedian(bbpProfile,5,align = "center", "extend"),5,align = "center","extend")
        
        bbpPresTemp<-approx(presProfile,bbpProfile,
                         1:max(presProfile))$x
        bbpTemp<-approx(presProfile,bbpProfile,
                            1:max(presProfile))$y
        bbpSmoothTemp<-approx(presProfile,bbpSmoothProfile,
                              1:max(presProfile))$y
        
        #Minimum values should be based on unsmoothed interpolated profile,
        #since minimum values are not affected by spikes
        #First identify minimum depth and value within top 10 meters
        bmin1_index <- floor(which.min(bbpSmoothTemp[1:10]))
        bmin1 = bbpSmoothTemp[bmin1_index]
        #Then identify the maximum beneath the top 10 meters (e.g., BBP_MAX),
        #based on the smoothed profile
        bmax1_index <- floor(which.max(bbpSmoothTemp[10:length(bbpSmoothTemp)])) + 10
        bmax1 = bbpSmoothTemp[bmax1_index]
        #Then identify the depths of minima beneath the maximum value, selecting the first
        #value as the BBP minimum
        bmin2_index <- rollapply(as.zoo(bbpSmoothTemp[bmax1_index:length(bbpSmoothTemp)]),
                  61, partial = FALSE, align='center',
                  function(x) which.min(x)==31)
        bmin2_index<-index(bmin2_index)[coredata(bmin2_index)][1] + bmax1_index
        
        #If no minima are identified, skip the profile
        if(is.na(bmin2_index)==TRUE){
          next
        }
      
        #Determine deep BBP minimum value based on the unsmoothed profile
        bmin2 = bbpSmoothTemp[bmin2_index]
        
        #Determine maximum value beneath the deep BBP minimum (e.g., OMZ maximum)
        #Again, use the smoothed profile
        bmax2_index <- which.max(bbpSmoothTemp[bmin2_index:length(bbpSmoothTemp)]) + bmin2_index
        bmax2 = bbpSmoothTemp[bmax2_index]
        #bmed = median(bbpTemp[1:10],na.rm=TRUE)

        temp$BBP_MAX[range] = bmax1
        temp$BBP_MAX_DEPTH[range] = bmax1_index
        temp$BBP_MESO_MAX[range] =  bmax2
        temp$BBP_MESO_MAX_DEPTH[range] = bmax2_index
        
        temp$DCM_BBP[range] = bbpSmoothTemp[cmax_index]
        temp$DCM_BBP_RATIO[range] = chlSmoothTemp[cmax_index]/bbpSmoothTemp[cmax_index]
        temp$BBP_MAX_MIN[range] = round(bmax1/bmin1,2)
        temp$BBP_MESO_MAX_MIN[range] = round(bmax2/bmin2,2)
        
        bbpProfilePlot<-ggplot()+
          scale_y_continuous(limits=c(0,plotMaxBBP))+
          coord_flip() + 
          scale_x_reverse(limits=c(750,0)) + 
          geom_line(aes(x=bbpPresTemp, y=bbpTemp),linetype = 2, alpha = 0.3) + 
          geom_line(aes(x=bbpPresTemp, y=bbpSmoothTemp)) + 
          geom_vline(aes(xintercept = bmax1_index),col="black") + 
          geom_vline(aes(xintercept = bmax2_index),col="black") + 
          geom_vline(aes(xintercept=bmin2_index),linetype=2) + 
          theme_bw() + 
          theme(plot.title = element_text(hjust=0.5),
                axis.title.y = element_blank()) +
          labs(y = bquote(b[bp]~"("*m^-1*")"))
          
          #Add green line 
          if(isDCM == TRUE &
             bmax1_index <= (cmax_index + 20) & 
             bmax1_index >= (cmax_index - 20) &
             round(bmax1/bmin1,2) >= 1.3){
              bbpProfilePlot <- bbpProfilePlot + 
                geom_vline(xintercept = bmax1_index,col ="#91A737")
          }
        
        if(round(bmax2/bmin2,2) >= 1.15){
          bbpProfilePlot <- bbpProfilePlot + 
            geom_vline(xintercept = bmax2_index,col ="#91A737")}
        }
      
      if(plot==TRUE){
        currentMonth <- month(currentDate)
        plotRow = plot_grid(chlProfilePlot,bbpProfilePlot,ncol=2)
          
        title <- ggdraw() + 
          draw_label(
              paste("Chl = ",round(cmax/cmed,2),"|||",
                    "BBP_MAX =",round(bmax1/bmin1,2),"|||",
                    "BBP_MESO =",round(bmax2/bmin2,2),"|||",
                    "i =",i,"j =",j),
              fontface = 'bold',
              x = 0,
              hjust = -.45
            ) +
            theme(
              # add margin on the left of the drawing canvas,
              # so title is aligned with left edge of first plot
              plot.margin = margin(0, 0, 0, 7))
          
          output<-plot_grid(
            title, plotRow,
            ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1))
          

          ggsave(paste("figures/parameterTestPlots/dcmTest/",
                       str_pad(currentMonth,width = 2, pad = 0),"_",
                       currentFloatID,"_",str_pad(j,width=3,pad=0),".png",sep=""),plot=output,
                 width = 10, height = 5,
                 dpi = 300)
        }
    }
    
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
    }
}

bathDepth<-function(fileList,replace=FALSE,
                    inputDirectory = "data/finalCSV_params/",
                    outputDirectory = "data/finalCSV_params/"){
  #browser()
  require(RNetCDF)
  bathymetry <- read.nc(open.nc("/Users/nicholasbock/Projects/argo/bathymetry/GEBCO_2020.nc"))
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep=""))
    if("bathDepth" %in% colnames(temp) & replace == FALSE){
      next
    }
    
    temp<-temp %>% 
      group_by(date) %>% 
      mutate(bathDepth = bathGet(lat=mean(latitude,na.rm=TRUE),
                                 lon=mean(longitude,na.rm=TRUE),
                                 bath = bathymetry))
    
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
                                
  }
  rm(bathymetry)
}

surfaceValues <- function(fileList,replace=FALSE,
                          inputDirectory = "data/finalCSV_params/",
                          outputDirectory = "data/finalCSV_params/"){
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep="")) %>% 
      mutate(PRES = as.numeric(PRES)) %>% 
      arrange(jDay,PRES)
    for(j in 1:length(unique(temp$date))){
      print(j)
      range = temp$date==unique(temp$date)[j]
      temp$SURF_TEMP[range] = mean(temp$TEMP[range][temp$PRES[range]<=5],na.rm=TRUE)
      temp$SURF_PSAL[range] = mean(temp$PSAL[range][temp$PRES[range]<=5],na.rm=TRUE)
      temp$SURF_SIG[range] = gsw_sigma0(SA = temp$SURF_PSAL[range][1],
                                        CT = temp$TEMP[range][1])
    }
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }
}

isolumes <- function(fileList,parFiles,
                     inputDirectory = "data/Chla_BBP_data/finalCSV/",
                     outputDirectory = "data/Chla_BBP_data/finalCSV/"){
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

oxycline <- function(fileList,plot=FALSE,
                     inputDirectory = "data/finalCSV/",
                     outputDirectory = "data/finalCSV/"){
  library(zoo)
  library(tidyverse)
  for(i in 1:length(fileList)){
    print(paste(i,"of",length(fileList)))
    temp<-read.csv(paste(inputDirectory,fileList[i],sep=""))
    currentFloat = temp$float[1]
    currentBasin = temp$basin[1]
    
    if("DOXY_DCM" %in% colnames(temp)){
      temp<-temp %>% 
        select(-DOXY_DCM,-DOXY_BBP_MAX)
    }
    
    temp$DOXY_MIN_DEPTH = NA
    temp$DOXY_MAX_DEPTH = NA
    temp$DOXY_BBP_MAX = NA
    temp$BBP_MAX_DOXY = NA
    temp$BBP_MESO_MAX_DOXY = NA
    temp$DCM_DOXY = NA
    temp$DOXY_MAX = NA
    temp$DOXY_MIN = NA
    temp$DOXY_SAT_MIN_DEPTH = NA
    temp$DOXY_SAT_MIN = NA
    #temp$DOXY_CLINE = NA
    #temp$DOXY_REMIN = NA
    if("DOXY" %in% colnames(temp)==FALSE){
      write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
                row.names = FALSE)
      next}
    for(j in 1:length(unique(temp$date))){
      print(j)
      #Identify current date
      currentDate = unique(temp$date)[j]
      #Identify rows in data in temp matching current date
      range = which(temp$date==currentDate)
      if(all(is.na(temp$DOXY[range]))==TRUE){
        next
      }
      
      doxyTemp <- na.fill(approx(temp$PRES[range],temp$DOXY[range],1:max(temp$PRES[range]))$y,"extend")
      presTemp <- approx(temp$PRES[range],temp$DOXY[range],1:max(temp$PRES[range]))$x
      if(all(is.na(temp$PSAL[range]))==FALSE & all(is.na(temp$TEMP[range]))==FALSE){
        tempTemp <- approx(temp$PRES[range],temp$TEMP[range],1:max(temp$PRES[range]))$y
        salTemp <- approx(temp$PRES[range],temp$PSAL[range],1:max(temp$PRES[range]))$y
        sat_data<-oxySat(doxyTemp,tempTemp,salTemp)
        satMinIndex = floor(which.min(sat_data$sat))
        temp$DOXY_SAT_MIN[range] = sat_data$sat[satMinIndex]
        temp$DOXY_SAT_MIN_DEPTH[range] = satMinIndex
      }
      
      #Identify maximum oxygen concentration (may be subsurface in regions w/ DCMs)
      doxyMinIndex = floor(which.min(doxyTemp))
      doxyMaxIndex = floor(which.max(doxyTemp))
      doxyMin = doxyTemp[doxyMinIndex]
      doxyMax = doxyTemp[doxyMaxIndex]

      #Subset data to only include depths greater than doxyMaxDepth
      # doxyTemp_subset<-doxyTemp[doxyMaxIndex:length(doxyTemp)]
      # presTemp_subset<-presTemp[doxyMaxIndex:length(doxyTemp)]
  
      #doxyThreshold = doxyMin + ((doxyMax - doxyMin)/2)
      #Identify depth where oxygen is for the first time is 10 um lower than max
      #doxyCline_mid_depth <- which.min(abs(doxyTemp-doxyThreshold))
      #doxyCline_start_depth <- floor(doxyCline_mid_depth*.8)
      #doxyCline_end_depth <- floor(doxyCline_mid_depth*1.25)
      #doxyCline_start_depth <- presTemp_subset[which((doxyTemp_subset - doxyMax + doxyThreshold) < 1)[1]]

      # doxyclineReg <- (summary(lm(doxyTemp[doxyCline_start_depth:doxyCline_end_depth]~
      #                                          presTemp[doxyCline_start_depth:doxyCline_end_depth])))
      #factor = factorList[which(rSq == max(rSq,na.rm=TRUE))]
      #factor = which(slope == median(slope,na.rm=TRUE))
      #doxyclineReg <- summary(lm(doxyTemp[doxyCline_start_depth:(doxyCline_start_depth)]~presTemp[doxyCline_start_depth:(doxyCline_start_depth)]))
      #doxyReminReg <- summary(lm(doxyTemp[doxyMinIndex:max(presTemp,na.rm=TRUE)]~presTemp[doxyMinIndex:max(presTemp,na.rm=TRUE)]))
      temp$DOXY_MIN_DEPTH[range] = doxyMinIndex
      temp$DOXY_MAX_DEPTH[range] = doxyMaxIndex
      temp$DOXY_MIN[range] = doxyMin
      temp$DOXY_MAX[range] = doxyMax
      temp$BBP_MAX_DOXY[range] = doxyTemp[temp$BBP_MAX_DEPTH[range][1]]
      temp$BBP_MESO_MAX_DOXY[range] = doxyTemp[temp$BBP_MESO_MAX_DEPTH[range][1]]
      temp$DCM_DOXY[range] = doxyTemp[median(temp$DCM_DEPTH[range])]
      #temp$DOXY_CLINE[range] = doxyReminReg$coefficients[1,2]
      #temp$DOXY_REMIN = doxyReminReg$coefficients[1,2]
      
      
      if(plot==TRUE){
        output<-ggplot() + 
          geom_line(aes(presTemp,doxyTemp),col="grey") + 
          geom_vline(xintercept = doxyMinIndex,col="grey",linetype = 2)+
          geom_vline(xintercept = doxyMaxIndex,col="grey",linetype = 2)+
          # geom_abline(aes(slope = -doxyclineReg$coefficients[2,1],
          #                 intercept = doxyclineReg$coefficients[1,1]),col = "red", linetype = 2) +
          # geom_abline(aes(slope = -doxyReminReg$coefficients[2,1],
          #                 intercept = doxyReminReg$coefficients[1,1]),col = "blue", linetype = 2) +
          scale_x_reverse() + 
          coord_flip() +
          theme_bw() + 
          labs(x = "PRES", y = "O2")
        
        if(all(is.na(temp$PSAL[range]))==FALSE & all(is.na(temp$TEMP[range]))==FALSE){
          output <- output + 
            geom_line(aes(presTemp,sat_data$sat),col="blue") +
            geom_vline(xintercept = satMinIndex,col="blue",linetype = 2)+ 
            labs(x = "PRES", y = "%Saturation (Blue); O2 (Grey)")
        }
        ggsave(plot = output,
               paste("figures/parameterTestPlots/oxyclineTest/",
               currentBasin,"_",currentFloat,"_",j,".png",sep=""))
      }
    }
    
    write.csv(temp, paste(outputDirectory,fileList[i],sep=""),
              row.names = FALSE)
  }
  
}

fileList<-list.files("data/finalCSV/",
                     full.names = F,pattern=".csv")
parFiles<-list.files("data/parData",pattern=".csv",
                      full.names = T)
# # #######RUN BEFORE PERFORMING PAR ESTIMATES#######
province(fileList)
mld(fileList,plot=FALSE)
####ChlSat takes about a day
#chlSat(fileList,replace=FALSE)
####canyonB takes about a day
#canyonB(fileList,replace=FALSE)
####intNewts takes about an hour
#intNewts(fileList)
###nclines takes about an hour or so
#nclines(fileList,plot=TRUE)
###DCM takes about an hour or so
#dcmDepth(fileList,plot=FALSE)
###Bath depth takes about 2 hours...there's probably a better way...
#bathDepth(fileList,replace=FALSE)
#surfaceValues(fileList)
###
#oxycline(fileList,plot=FALSE)

###########RUN AFTER PERFORMING PAR ESTIMATES###########
isolumes(fileList,parFiles)
bbpInt(fileList)


