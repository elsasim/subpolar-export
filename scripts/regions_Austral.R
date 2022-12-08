

floats <- read_csv("data/weightings_clusters_Austral.csv") 
floats <- as.vector(floats$floatID) %>%
  str_remove("_[0-9]") %>%
  unique()

time_series <- read_csv("data/weightings_clusters_Austral.csv") 
time_series <- as.vector(time_series$floatID)

for (i in 1:length(floats)){
  print(paste("File",i,"of",length(floats),sep = " "))
  file <- floats[i]
  import <- read_csv(paste("data/",file,".csv",sep=""), show_col_types = F) %>%
    select(floatID,latitude,longitude,province) %>%
    distinct() %>%
    filter(floatID %in% time_series)
  
  if (exists("final")){
    final <- bind_rows(final, import)
  }
  
  if (!exists("final")){
    final <- import
  }
  
}

final_ASZ_SIZ <- filter(final,province == "ASZ_SIZ")
final_PFZ <- filter(final,province == "PFZ")
final_ASZ_SIZ_PFZ <- filter(final,province == "ASZ_SIZ" | province == "PFZ")
final_SAZ <- filter(final,province == "SAZ")
final_STZ <- filter(final,province == "STZ")

library(rworldmap)
library(factoextra)
library(FactoMineR)
library(sf)
worldMap<-getMap(resolution="coarse")

clusterMap <- ggplot()+
  geom_polygon(data = worldMap, aes(x = long, y = lat, group=group),
               fill = "grey", col = "black") + 
  coord_map(ylim = c(-70, -25), xlim = c(-30, 60)) +
  geom_point(data = final,
             aes(longitude, latitude, fill = floatID), col = "black", shape = 21, alpha = 0.8, size = 5)+
  # geom_hline(aes(yintercept = 40),linetype = 2) + 
  # geom_hline(aes(yintercept = -40),linetype =2) + 
  guides(fill=guide_legend(byrow = TRUE)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(face="bold", size = 7),
        plot.title = element_text(hjust = 0.5),
        legend.box.spacing = unit(0.1, "pt")) +
  scale_fill_manual(values = c("#FF6666","#CC0000","#FFCC00"),
                    labels = c("5904469_2 (2)", "5904469_3 (1)", "5904469_4 (3)"))

ggsave(plot=clusterMap,
       "figures/non_cross_regions_map.png",
       width = 8, height = 5, dpi = 300)


final_2 <- group_by(final,floatID) %>%
  filter(length(unique(province)) != 1) %>%
  summarise(latitude = mean(latitude), longitude = mean(longitude), province = unique(province))




all_floats <- read_csv("data/weightings_clusters_Austral.csv") %>%
  filter(floatID %in% data$floatID ==F)

data_crossed <- final_2 %>%
  group_by(floatID)%>%
  summarise(latitude = first(latitude), longitude = first(longitude)) %>%
  mutate(cluster = as.character(all_floats$cluster))

data_non_crossed <- final %>%
  filter(floatID %in% data_crossed$floatID ==F) %>%
  group_by(floatID) %>%
  summarise(latitude = first(latitude), longitude = first(longitude)) %>%
  mutate(cluster = as.character(all_floats$cluster))

clusters <- read_csv("data/weightings_clusters_HL_N.csv") %>%
  filter(cluster == 2)



