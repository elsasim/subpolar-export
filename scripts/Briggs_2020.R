floats <- clusters$floatID

floats <- unique(str_extract(floats,"([:digit:]{7})"))

files <- list.files("data/Chla_BBP_data/files_POC/",pattern=".csv",full.names = F)
files <- str_extract(files,"([:digit:]{7})")

for (i in 1:length(files)){
  print(paste("File",i,"of",length(files)))
  import <- read_csv(files[i],show_col_types = F)
  
  currentFloat <- import$float[1]
  
  if(currentFloat %in% floats ==F){
    file.remove(paste("data/Chla_BBP_data/files_POC/",currentFloat,".csv",sep = ""))
  }
  
}
##############################################################################
##############################################################################
##############################################################################

# determination du bruit du capteur pour chaque flotteur

  
# appliquer roll filter min puis max pour avoir bbl + noise
    NOISE <- import %>% 
      select("PRES","date","cycleNumber","CHLA") %>%
      filter(is.na(CHLA)==FALSE) %>%
      group_by(cycleNumber) %>%
      mutate(nRows = length(cycleNumber)) %>%
      filter(nRows>1) %>%
      mutate(filtered = rollapply(CHLA, 11, min, fill = "extend")) %>%
      mutate(filtered = rollmax(filtered, 11, fill = "extend")) %>%
      filter(PRES>300) %>%
      ungroup() %>%
      mutate(bbl_noise = CHLA - filtered) %>%
      mutate(mediane = median(bbl_noise)) %>%
      
      # filter((bbl_noise<0)==F) %>%
      mutate (bin = NA)
    
    # ranger par tranche de 50m sous 300m
    for(j in 1:length(NOISE$CHLA)){
      for(i in 1:14){
        if(NOISE$PRES[j]<(300+i*50) & NOISE$PRES[j]>=(300+(i-1)*50)){
          NOISE$bin[j] <- i
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

    import <- import %>%
      mutate(noise=noise)
  
##############################################################################
##############################################################################
##############################################################################

# post traitement de la chla
if(any(import$CHLA<0, na.rm=T)){

    for(j in 1:length(unique(import$cycleNumber))){
      cycle <- import %>%
        filter(cycleNumber == j)
      if(any(cycle$CHLA<0, na.rm=T)){
        CHLA <- cycle$CHLA
        CHLA_1 <- CHLA + abs(min(CHLA, na.rm = T)) + noise #si valeurs negatives alors mettre la valeur min à la valeur du bruit
        import$CHLA[which(import$cycleNumber==j)] <- CHLA_1
        print(paste("cycle :",j))
      }
    }
}
