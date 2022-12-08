require(tidyverse)
require(lubridate)

fileList<-list.files("/Users/nicholasbock/Projects/argo/data/finalCSV_old",
                     full.names = T)
#dcmDepth(fileList,plot=TRUE)

Sys.setenv(tz = "UTC")
fileList<-list.files("/Users/nicholasbock/Projects/argo/data/finalCSV_old/",
                     full.names = T)

clusterList<-read.csv("/Users/nicholasbock/Projects/argo/annualSVD/dataCompilations/finalClusters.csv") %>% 
  select(floatID,cluster=sixClusters)

outputFile <- read.csv(fileList[7]) %>% 
  mutate(month = month(date),
         jDay = yday(date),
         month_adj = (((jDay-1)*12)/(max(jDay)-1))+1,
         month = factor(month)) %>% 
  merge(clusterList,by="floatID",all.x=TRUE) %>% 
  filter(is.na(cluster)==FALSE)

#group by parameters that are unique to each profile
outputFile<-outputFile %>% 
  select(cluster,float,timeSeries,floatID,basin,province,date,date_orig,jDay,
         latitude,longitude,bathDepth,
         MLD_003, BVF,
         CHL_SAT, SURF_TEMP, SURF_PSAL, SURF_SIG,
         DCM, DCM_DEPTH,DCM_MED_MAX, DCM_DIFF, DCM_BBP,DCM_BBP_RATIO,
         BBP_MAX,BBP_MAX_DEPTH,BBP_MED_MAX,
         AVG_MLD_NO3, AVG_MLD_CAN_NO3, AVG_MLD_CAN_PO4, AVG_MLD_CAN_SiO4,
         NCLINE_DEPTH, NCLINE_SLOPE, NCLINE_R2, 
         CAN_NCLINE_DEPTH, CAN_NCLINE_SLOPE,CAN_NCLINE_R2,
         CAN_PCLINE_DEPTH, CAN_PCLINE_SLOPE,CAN_PCLINE_R2,
         CAN_SCLINE_DEPTH, CAN_SCLINE_SLOPE,CAN_SCLINE_R2,
         ISO_10, ISO_1, ISO_01,IPAR20,PAR_FLAG,
         NCLINE_PAR, CAN_NCLINE_PAR,CAN_PCLINE_PAR, CAN_SCLINE_PAR,MLD_MED_PAR,
         AVG_BBP_MESO_0_4,INT_BBP_MESO_0_4,
         AVG_BBP_MESO_0_1,AVG_BBP_MESO_1_2,AVG_BBP_MESO_2_3,AVG_BBP_MESO_3_4,
         INT_BBP_MESO_0_1,INT_BBP_MESO_1_2,INT_BBP_MESO_2_3,INT_BBP_MESO_3_4,
         INT_BBP_PZ, INT_BBP_EZ, INT_CHL_PZ,INT_CHL_EZ,
         INT_BBP_MESO_YR,INT_BBP_PZ_YR,
         AVG_BBP_PZ_YR,
         AVG_BBP_MESO_SLOPE,AVG_BBP_MESO_SLOPE_R) %>% 
    group_by(cluster,float,timeSeries,floatID,basin,province,date,date_orig,jDay,latitude,longitude,bathDepth) %>% 
    summarize_all(mean,na.rm=TRUE)

for(i in 8:length(fileList)){
  print(i)
  temp <- read.csv(fileList[i])  %>% 
    mutate(month = month(date),
           jDay = yday(date),
           month_adj = (((jDay-1)*12)/(max(jDay)-1))+1,
           month = factor(month)) %>% 
    merge(clusterList,by="floatID",all.x=TRUE) %>% 
    filter(is.na(cluster)==FALSE)
  
  if(nrow(temp)==0){
    print("blegh")
    next
  }


  temp<-temp %>% 
    #group by parameters that are unique to each profile
    select(cluster,float,timeSeries,floatID,basin,province,date,date_orig,jDay,
           latitude,longitude,bathDepth,
           MLD_003, BVF,
           CHL_SAT, SURF_TEMP, SURF_PSAL, SURF_SIG,
           DCM, DCM_DEPTH,DCM_MED_MAX, DCM_DIFF, DCM_BBP,DCM_BBP_RATIO,
           BBP_MAX,BBP_MAX_DEPTH,BBP_MED_MAX,
           AVG_MLD_NO3, AVG_MLD_CAN_NO3, AVG_MLD_CAN_PO4, AVG_MLD_CAN_SiO4,
           NCLINE_DEPTH, NCLINE_SLOPE, NCLINE_R2, 
           CAN_NCLINE_DEPTH, CAN_NCLINE_SLOPE,CAN_NCLINE_R2,
           CAN_PCLINE_DEPTH, CAN_PCLINE_SLOPE,CAN_PCLINE_R2,
           CAN_SCLINE_DEPTH, CAN_SCLINE_SLOPE,CAN_SCLINE_R2,
           ISO_10, ISO_1, ISO_01,IPAR20,PAR_FLAG,
           NCLINE_PAR, CAN_NCLINE_PAR,CAN_PCLINE_PAR, CAN_SCLINE_PAR,MLD_MED_PAR,
           AVG_BBP_MESO_0_4,INT_BBP_MESO_0_4,
           AVG_BBP_MESO_0_1,AVG_BBP_MESO_1_2,AVG_BBP_MESO_2_3,AVG_BBP_MESO_3_4,
           INT_BBP_MESO_0_1,INT_BBP_MESO_1_2,INT_BBP_MESO_2_3,INT_BBP_MESO_3_4,
           INT_BBP_PZ, INT_BBP_EZ, INT_CHL_PZ,INT_CHL_EZ,
           INT_BBP_MESO_YR,INT_BBP_PZ_YR,
           AVG_BBP_PZ_YR,
           AVG_BBP_MESO_SLOPE,AVG_BBP_MESO_SLOPE_R)  %>% 
    group_by(cluster,float,timeSeries,floatID,basin,province,date,date_orig,jDay,latitude,longitude,bathDepth) %>% 
  summarize_all(mean,na.rm=TRUE)
  
  outputFile<-rbind(outputFile,temp)
}

outputFile<-outputFile %>% 
  mutate(month_adj = (((jDay-1)*12)/(max(jDay)-1))+1)

outputFile$qtr = NA
outputFile$qtr[month(outputFile$date)%in%c(12,1,2)] = 1
outputFile$qtr[month(outputFile$date)%in%c(3,4,5)] = 2
outputFile$qtr[month(outputFile$date)%in%c(6,7,8)] = 3
outputFile$qtr[month(outputFile$date)%in%c(9,10,11)]= 4

outputFile$qtr[outputFile$cluster=="AR" & outputFile$month_adj >= 12] = 1
outputFile$qtr[outputFile$cluster=="AR" & outputFile$month_adj < 4] = 1
outputFile$qtr[outputFile$cluster=="AR" & outputFile$month_adj >= 4 & outputFile$month_adj < 7] = 2
outputFile$qtr[outputFile$cluster=="AR" & outputFile$month_adj >= 7 & outputFile$month_adj < 9.5] = 3
outputFile$qtr[outputFile$cluster=="AR" & outputFile$month_adj >= 9.5 & outputFile$month_adj < 12] = 4
outputFile$qtr <- factor(outputFile$qtr,levels=c(1,2,3,4))

write_csv(outputFile,"/Users/nicholasbock/Projects/argo/annualSVD/dataCompilations/allParams_0701.csv")
