

POC_sd


layers_names <- c("0_800","100_800", "200_800", "300_800", "400_800", "500_800", "600_800", "700_800")
incertitudes <- data.frame(cluster = c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8)), layer = rep(layers_names,5), slope_min = rep(NA,40), slope_max = rep(NA,40))



for (i in 1:8){

min_y <- get(paste("min_",layers_names[i],"_y", sep=""))
max_y <- get(paste("max_",layers_names[i],"_y", sep=""))
min_x <- get(paste("min_",layers_names[i],"_x", sep=""))
max_x <- get(paste("max_",layers_names[i],"_x", sep=""))
sd_tot <- get(paste("sd_tot_",layers_names[i], sep=""))[[5]]

if(min_x>max_x){
    # pente minimum
    new_min_y <- min_y + sd_tot[min_x]
    new_max_y <- max_y - sd_tot[max_x]
    slope_min <- (new_max_y - new_min_y) / (max_x - min_x) 
    
    # pente maximum
    new_min_y <- min_y - sd_tot[min_x]
    new_max_y <- max_y + sd_tot[max_x]
    slope_max <- (new_max_y - new_min_y) / (max_x - min_x) 
}

if(min_x<max_x){
    # pente minimum
    new_min_y <- min_y + sd_tot[min_x]
    new_max_y <- max_y - sd_tot[max_x]
    slope_min <- (new_min_y - new_max_y) / ((365-max_x) + min_x)

    # pente maximum
    new_min_y <- min_y - sd_tot[min_x]
    new_max_y <- max_y + sd_tot[max_x]
    slope_max <- (new_min_y - new_max_y) / ((365-max_x) + min_x)
}


incertitudes$slope_min[i+32] <- -(1/1.45)*slope_min*10^-3 *12
incertitudes$slope_max[i+32] <- -(1/1.45)*slope_max*10^-3 *12
# incertitudes$slope_min[i+32] <- slope_min *10^3
# incertitudes$slope_max[i+32] <- slope_max *10^3

}
write_csv(incertitudes,"data/all/incertitudes_O2.csv")

