setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig1_nj_structure_he_plastome/Fig1PCA/")

library(mapdata)
library(usmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(sp)
library(ggspatial)
library(ggmap)
library(ggsn)
library(dplyr)
library(maps)
library(maptools)
library(ggplot2)
library(ggrepel)
library(grid)
library(cowplot)

##################################################################
############## Caribbean region ##################################
##################################################################

register_google("AIzaSyAFMCC4ztGaZzBKvTMo4ByO5pBZVPhR5tQ") #for the ggmap package

gpsmk <- read.csv("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig0_map/CottonSamplingMap5_ut8.csv", encoding="UTF-8")


lon_range <- c(-83.4, -79.5)  # Longitude (xlim)
lat_range <- c(24.05, 28.2)     # Latitude (ylim)

################################

gpsmk2 <- gpsmk %>%
  filter(Source == "Field") %>%
  filter(Population %in% c("Yucatan", "Florida")) %>%
  filter(Group != "Landrace") %>%
  mutate(maploc = paste0(
    stringr::str_extract(PopInfor, "^.{2}"),  stringr::str_extract(Population2, "^.{2}") 
    #" (",
    #PopInfor, "_", Population2, 
    #")"
  )) %>%
  add_count(maploc, name = 'id_occurence')  %>%
  mutate(maploc2 = paste0( maploc, " (n = ", id_occurence, ")"))

################################
gpsmk3 <- gpsmk %>%
  filter(Source == "Field") %>%
  filter(Population == "Yucatan")  %>%
  filter(Group != "Landrace") %>%
  mutate(maploc = paste0(
    stringr::str_extract(PopInfor, "^.{2}"),  stringr::str_extract(Population2, "^.{2}") 
    #" (",
    #PopInfor, "_", Population2, 
    #")"
    )) %>%
  add_count(maploc, name = 'id_occurence')  %>%
  mutate(maploc2 = paste0( maploc, " (n = ", id_occurence, ")"))

################################

gpsmk4 <- gpsmk %>%
  filter(Source == "Field") %>%
  filter(Region == "Florida")  %>%
  filter(Group != "Landrace") %>%
  add_count(Group, name = 'id_occurence')  %>%
  mutate(maploc = paste0( Group, " (n = ", id_occurence, ")")) %>%
  mutate(maploc = gsub("MK", "FL_MK", maploc))   %>%
  mutate(maploc = gsub("AD1_", "", maploc)) %>%
  mutate(maploc = gsub("FL_", "", maploc))


##################################################################
##################################################################
##################################################################
# Get the map from Google
map1 <- get_map(c(left = -98, bottom = 4, right = -55, top = 31), 
               zoom = 4, 
               source = "google", 
               maptype = "satellite")

# Plot the map
bigmap <- ggmap(map1) +
  xlim(lon_range) +
  ylim(lat_range) +
  coord_sf(xlim = c(-98, -55), ylim = c(4, 31), expand = FALSE, crs = "WGS84")+
  xlab("Longitude") + ylab("Latitude") +
  geom_point(data = gpsmk2, mapping = aes(x = long, y = lat-2,  shape = Source), color = "yellow", size =4, show.legend = F) +
  annotate(geom = "text", x = -84, y = 25, label = "Mound Key", fontface = "bold", color = "white", size = 4) +
  annotate(geom = "text", x = -65, y = 17, label = "Puerto Rico", fontface = "bold", color = "white", size = 4) +
  annotate(geom = "text", x = -60, y = 15, label = "Guadeloupe", fontface = "bold", color = "white", size = 4) +
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+
  scalebar( x.min = -94, x.max = -90, y.min = 6, y.max = 8,
            box.fill = c("yellow", "white"), st.color = "white",
            dist = 200, dist_unit = "km", st.size =3, height =0.3, st.dist = 0.2,
            transform = TRUE, model = "WGS84") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.15, "in"),
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                line_col = 'white')) +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(family = 'serif'),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.1),
        panel.grid.major = element_line(color = "grey90", linetype = "dashed",  size = 0.2))




####################################################################
####################################################################
####################################################################

world <- ne_countries(scale = "medium", type = "countries", returnclass = "sf")
bigmap2 <- ggplot(data = world) +
  geom_sf(fill = "#78909C", color = "white", show.legend = F) +
  coord_sf(xlim = c(-98, -55), ylim = c(6, 31), expand = FALSE)+
  #xlab("Longitude") + ylab("Latitude") +
 # annotate(geom = "text", x = -65, y = 19, label = "Puerto Rico", fontface = "bold", color = "black", size = 3) +
 # annotate(geom = "text", x = -60, y = 17, label = "Guadeloupe", fontface = "bold", color = "black", size = 3) +
  
  scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01))+
  
theme(rect = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent",  colour = NA_character_),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      panel.background = element_rect(fill = "white", colour = NA_character_), # necessary to avoid drawing panel outline
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.title=element_blank())
bigmap2
















##################################################################
##################################################################
##################################################################


# Get the map from Google
map <- get_map(location = c(lon = mean(lon_range), lat = mean(lat_range)), 
               zoom = 7, 
               source = "google", 
               maptype = "satellite")

ggmap(map)

FL_samplingplot <- ggmap(map) +
  coord_sf(xlim = c(-83.4, -79.5), ylim = c(24.05, 28.2), expand = FALSE, crs = "WGS84") +
  geom_point(data = gpsmk4, mapping = aes(x = long, y = lat-0.04), color = "yellow", size = 3, show.legend = F) +
  geom_text_repel(data = gpsmk4[!duplicated(gpsmk4[,"maploc"]),], aes(x = long, y = lat -0.04, label=maploc), 
                  min.segment.length = unit(0, 'lines'),  color = "white", 
                  max.overlaps=Inf, size = 5, show.legend = F ) +
  scalebar( x.min = -83, x.max = -82, y.min = 24.3, y.max = 25,
            box.fill = c("yellow", "white"), st.color = "white",
            dist = 50, dist_unit = "km", st.size =3, height =0.14, st.dist = .12,
            transform = TRUE, model = "WGS84") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                line_col = 'white')) +
  theme_bw() +
  theme(rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_rect(fill = "white", colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title=element_blank())

FL_samplingplot

##################################################################
##################################################################
##################################################################

lon_range <- c(-91, -86.5)  # Longitude (xlim)
lat_range <- c(20, 22)     # Latitude (ylim)

# Get the map from Google
map2 <- get_map(location = c(lon = mean(lon_range), lat = mean(lat_range)), 
               zoom = 7, 
               source = "google", 
               maptype = "satellite")

ggmap(map2)

YUC_samplingplot <- ggmap(map2) +
  coord_sf(xlim = c(-91, -86.5), ylim = c(20, 22.5),  expand = FALSE, crs = "WGS84") +
  geom_point(data = gpsmk3, mapping = aes(x = long, y = lat-0.04), color = "yellow", size = 3, show.legend = F) +
  geom_text_repel(data = gpsmk3[!duplicated(gpsmk3[,"maploc2"]),], aes(x = long, y = lat -0.04, label=maploc2), 
                  min.segment.length = unit(0.3, 'lines'),  color = "white", 
                  max.overlaps=Inf, size = 5, show.legend = F ) +
  scalebar( x.min = -90.5, x.max = -90, y.min = 20.2, y.max = 20.4,
            box.fill = c("yellow", "white"), st.color = "white",
            dist = 25, dist_unit = "km", st.size =3, height =0.3, st.dist = 0.2,
            transform = TRUE, model = "WGS84") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering(text_col = 'white',
                                                                line_col = 'white')) +
  theme(rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        panel.background = element_rect(fill = "white", colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title=element_blank())

  
YUC_samplingplot


YUCinsert <- ggdraw() +
  draw_plot(YUC_samplingplot, 0, 0, 1, 1) +  # Main plot takes full area
  draw_plot(bigmap2, 0.7, 0.6, 0.3, 0.3)  # Inset plot in top-right corner

##################################################################
##################################################################
##################################################################
mapfinalplot <- ggdraw() +
  draw_plot(YUCinsert , x = 0, y = 0, width = 0.6, height = 1) +
  draw_plot(FL_samplingplot, x = 0.6, y = 0, width = 0.4, height = 1) 

#  draw_plot(bigmap, x = 0, y = 0.5, width = 0.6, height = 0.5) +
#  draw_plot(PR_samplingplot , x = 0.6, y = 0.75, width = 0.4, height = 0.25) +
#  draw_plot(GD_samplingplot, x = 0.6, y = 0.5, width = 0.4, height = 0.25) +
  
#  draw_plot(pca_plot, x = 0, y = 0, width = 0.6, height = 0.5) +
#  draw_plot(GDnjtree, x = 0.6, y = 0, width = 0.4, height = 0.5) +
#  draw_plot_label(label = c("b)", "c)"), size = 20,  fontface = "bold",
#                  x = c(0, 0.5), y = c(1, 1))

pdf("FL_YUC_map.pdf", width = 16, height = 8)
mapfinalplot
dev.off()

