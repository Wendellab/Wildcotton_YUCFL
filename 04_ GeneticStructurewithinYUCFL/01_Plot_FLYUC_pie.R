setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig5_PopGeneStructure_YUCFL/FigS678_nj_structure_relatedness")
library(geodata)
library(terra)  # For better raster handling
library(geodata)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(ggplot2)
library(dplyr)
library(LEA)
library(RColorBrewer)
library(stringr)
library(ggrepel)
library(scales)
library(ggpubr)

################################################################################################
gpsmk <- read.csv("supplementarytable.csv", encoding="UTF-8")
world <- ne_countries(scale = "large", returnclass = "sf")
world <- st_make_valid(world)

################################################################################################
################################################################################################


load("YUCn188_LEAstruture_n65.RData")

colplot <- plot_list[[11]]$data %>%
  mutate(individual = gsub("AD1_*", "", individual)) 
colnames(colplot)[4] = 'pop_bar'

YUC_pie <- colplot %>%
  pivot_wider(id_cols = c(individual, accession), names_from = pop_bar, values_from = q) %>%
  left_join(gpsmk %>%
              rename(individual = Individual) %>%
              mutate(individual = gsub("^AD1_", "", individual)) %>%
              select(individual, lat, long),
            by = "individual") %>%
  filter(!grepl("Cultivar|LR1|LR2", accession)) %>%
  mutate(accession = str_extract(individual, "^[^_]+_[^_]+")) %>%
  group_by(accession) %>%
  mutate(samplesnumber = n()) %>%
  group_by(accession, samplesnumber, lat, long) %>%
  summarise(across(P1:P11, mean, na.rm = TRUE), .groups = "drop")

#####################
# PCA-consistent colors for accessions
YUCcbbPalette <- colorRampPalette(brewer.pal(12, "Paired"))(8)
YUCcategories <- c("YUC_CeCo","YUC_CeDo","YUC_CeDr","YUC_CePr",
                   "YUC_RiCa","YUC_RiCh","YUC_SiPa","YUC_SiPr")
YUCcolormapping <- setNames(YUCcbbPalette, YUCcategories)

# Major component per accession + PCA color
YUC_major_component_per_accession <- YUC_pie %>%
  mutate(
    major_component = apply(select(., P1:P11), 1, function(x) names(x)[which.max(x)]),
    color = YUCcolormapping[accession]
  ) %>%
  select(accession, major_component, color)

# Assign colors to all P1–P11 components
YUC_all_components <- paste0("P", 1:11)
YUC_component_color_mapping <- setNames(RColorBrewer::brewer.pal(12, "Paired")[1:11], YUC_all_components)
YUC_component_color_mapping[YUC_major_component_per_accession$major_component] <- YUC_major_component_per_accession$color

YUC_component_color_mapping

######################

# Then crop
lon_range <- c(-90.7, -86.9)
lat_range <- c(19.5, 23)
YUCworld_crop <- st_crop(world,
                      xmin = lon_range[1], xmax = lon_range[2],
                      ymin = lat_range[1], ymax = lat_range[2])

YUC_pie_off <- YUC_pie %>%
  mutate(
    lat_plot  = lat,
    long_plot = long,
    lat_plot  = ifelse(accession == "YUC_CeCo", lat + 0.35, lat_plot),
    lat_plot  = ifelse(accession == "YUC_CeDr", lat - 0.35, lat_plot),
    
    lat_plot  = ifelse(accession == "YUC_CeDo", lat - 0.35, lat_plot),
    lat_plot  = ifelse(accession == "YUC_CePr", lat + 0.35, lat_plot),
    
    lat_plot  = ifelse(accession == "YUC_RiCa", lat + 0.35, lat_plot),
    lat_plot  = ifelse(accession == "YUC_RiCh", lat - 0.35, lat_plot),
    lat_plot  = ifelse(accession == "YUC_SiPa", lat + 0.35, lat_plot),
    lat_plot  = ifelse(accession == "YUC_SiPr", lat - 0.35, lat_plot))



# get country outlines

pie_cols <- paste0("P", 1:11)

YUCpie <- ggplot() +
  geom_sf(data = YUCworld_crop, fill = "gray95", color = "gray60", linewidth = 0.4) +
  geom_scatterpie(data = YUC_pie_off, aes(x = long_plot, y = lat_plot, r = sqrt(samplesnumber) * 0.04),
                  cols = pie_cols, color = NA) +
  geom_point(data = YUC_pie_off, aes(x = long, y = lat), color = "white", size = 3, shape = 24, fill = "black")+
  
  geom_text_repel(data = YUC_pie_off, aes(x = long, y = lat, label = gsub("^YUC_", "", accession)), 
                  max.overlaps = Inf,          # ← force all labels to show
                  size = 3, min.segment.length = 0, box.padding = 0.2, point.padding = 0.3) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", fill = "Ancestry") +
  scale_fill_manual(values = YUC_component_color_mapping, guide = guide_legend(ncol = 2))+
  
  scale_x_continuous(breaks = lon_range, labels = label_number(suffix = "°")) +
  scale_y_continuous(breaks = lat_range, labels = label_number(suffix = "°")) +
  theme_classic2() +
  theme(
    panel.background = element_blank(),    
    plot.background = element_blank(),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    #axis.title.x = element_text(size = 10, color = "black"),
    #axis.title.y = element_text(size = 10, color = "black", angle = 90),
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #panel.border = element_rect(color = "black", fill = NA),
    # Position legend inside plot
    legend.position = c(0.77, 0.3),  # Top-left corner
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA, color = "black", linewidth = 0.3),
    legend.title = element_text(size = 8, margin = margin(b = 2)),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.05, "cm"),  # Less vertical space
    legend.spacing.x = unit(0.1, "cm"),   # Space between columns
    legend.margin = margin(3, 3, 3, 3))

YUCpie  

###################################################################################
###################################################################################

load("FLn196_LEAstruture.RData")

colplot <- plot_list[[15]]$data %>%
  mutate(individual = gsub("AD1_*", "", individual)) 
colnames(colplot)[4] = 'pop_bar'

FL_pie <- colplot %>%
  pivot_wider(id_cols = c(individual, accession), names_from = pop_bar, values_from = q) %>%
  left_join(gpsmk %>%
              rename(individual = Individual) %>%
              mutate(individual = gsub("MK_Pop", "MK_Site", individual))%>%
              
              mutate(individual = gsub("^AD1_", "", individual)) %>%
              
              select(individual, lat, long),
            by = "individual") %>%
  filter(!grepl("Cultivar|LR1|LR2", accession)) %>%
  mutate(accession = str_extract(individual, "^[^_]+_[^_]+"),
         accession = case_when(
           accession %in% c("FL_OPB4", "FL_OPB5") ~ "FL_OPB4_OPB5",
           TRUE ~ accession)) %>%
  group_by(accession) %>%
  mutate(samplesnumber = n()) %>%
  group_by(accession, samplesnumber, lat, long) %>%
  summarise(across(P1:P15, mean, na.rm = TRUE), .groups = "drop")

############################################

# PCA-consistent colors for FL accessions
FLcbbPalette <- colorRampPalette(brewer.pal(12, "Paired"))(15)
FLcategories <- c("FL_CPH","FL_CPT","FL_FMY","FL_NP","FL_OPB2",
                  "FL_OPB4","FL_OPB5","FL_PK","FL_RBD","FL_RBT",
                  "FL_RNRB","FL_SR","FL_TC","FL_VKPA","FL_MK")
FLcolormapping <- setNames(FLcbbPalette, FLcategories)

# Major component per accession + assign PCA color
FL_major_component_per_accession <- FL_pie %>%
  mutate(
    major_component = apply(select(., P1:P15), 1, function(x) names(x)[which.max(x)]),
    color = FLcolormapping[accession]
  ) %>%
  select(accession, major_component, color)

# Assign colors to all P1–P15 components
FL_all_components <- paste0("P", 1:15)
FL_component_color_mapping <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(15), FL_all_components)
FL_component_color_mapping[FL_major_component_per_accession$major_component] <- FL_major_component_per_accession$color
FL_component_color_mapping["P1"] <- "#EC9A91"
FL_component_color_mapping["P15"] <- "#72ADD1"


FL_component_color_mapping


############################################

lon_range <- c(-83.6, -79.5)  # Longitude (xlim)
lat_range <- c(24.4, 28.2)     # Latitude (ylim)
FLworld_crop <- st_crop(world,
                      xmin = lon_range[1], xmax = lon_range[2],
                      ymin = lat_range[1], ymax = lat_range[2])

unique(FL_pie$accession)

FL_pie_off <- FL_pie %>%
  mutate(
    lat_plot  = lat,
    long_plot = long,
    
    # Earlier offsets
    long_plot = ifelse(accession == "FL_TC",  long + 0.3, long_plot),
    long_plot = ifelse(accession == "FL_NP",  long - 0.3, long_plot),
    long_plot = ifelse(accession == "FL_MK",  long - 0.3, long_plot),
    
    long_plot = ifelse(accession == "FL_CPT", long - 0.3, long_plot),
    long_plot = ifelse(accession == "FL_FMY", long - 0.4, long_plot),
    lat_plot  = ifelse(accession == "FL_FMY", lat  + 0.4, lat_plot),
    
    lat_plot  = case_when(
      accession == "FL_RBD"  ~ lat + 0.55,
      accession == "FL_RNRB" ~ lat + 0.4,
      accession == "FL_RBT"  ~ lat + 0.2,
      TRUE ~ lat_plot ),
    long_plot = case_when(
      accession %in% c("FL_RBD", "FL_RNRB", "FL_RBT") ~ long - 0.05,  # same x for all three
      TRUE ~ long_plot),
    
    long_plot = ifelse(accession == "FL_CPH", long - 0.3, long_plot),
    lat_plot  = ifelse(accession == "FL_VKPA", lat - 0.2, lat_plot),
    
    # ---- LAST FOUR SITES: one column with gradual x decrease ----
    lat_plot = case_when(
      accession == "FL_SR"        ~ lat + 0.2,
      accession == "FL_OPB2"      ~ lat + 0,
      accession == "FL_OPB4_OPB5" ~ lat - 0.30,
      accession == "FL_PK"        ~ lat - 0.45,
      TRUE ~ lat_plot),
    long_plot = case_when(
      accession == "FL_SR"        ~ max(long) + 0.4,
      accession == "FL_OPB2"      ~ max(long) + 0.2,   # slightly smaller x
      accession == "FL_OPB4_OPB5" ~ max(long) + 0,
      accession == "FL_PK"        ~ max(long) - 0.2,
      TRUE ~ long_plot))




levels(as.factor(pie_cols))


pie_cols <- paste0("P", 1:15)

FLCpie <- ggplot() +
  geom_sf(data = FLworld_crop, fill = "gray95", color = "gray60", linewidth = 0.4) +
  geom_scatterpie(data = FL_pie_off, aes(x = long_plot, y = lat_plot, r = sqrt(samplesnumber) * 0.04),
                  cols = pie_cols, color = NA) +
  geom_point(data = FL_pie_off, aes(x = long, y = lat), color = "white", size = 3, shape = 24, fill = "black")+
  geom_text_repel(data = FL_pie_off, aes(x = long, y = lat, label = gsub("^FL_", "", accession)), 
                  max.overlaps = Inf,          # ← force all labels to show
                  size = 3, min.segment.length = 0, box.padding = 0.2, point.padding = 0.3) +
  
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", fill = "Ancestry") +
  scale_fill_manual(values = FL_component_color_mapping, guide = guide_legend(ncol = 3))+
  
  scale_x_continuous(breaks = lon_range, labels = label_number(suffix = "°")) +
  scale_y_continuous(breaks = lat_range, labels = label_number(suffix = "°")) +
  theme_classic2() +
  theme(
    panel.background = element_blank(),    
    plot.background = element_blank(),
    #axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    #axis.title.x = element_text(size = 10, color = "black"),
    #axis.title.y = element_text(size = 10, color = "black", angle = 90),
    #axis.text.x = element_blank(),
    #axis.text.y = element_blank(),
    #panel.border = element_rect(color = "black", fill = NA),
    # Position legend inside plot
    legend.position = c(0.67, 0.95),  # Top-left corner
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA, color = "black", linewidth = 0.3),
    legend.title = element_text(size = 8, margin = margin(b = 2)),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.05, "cm"),  # Less vertical space
    legend.spacing.x = unit(0.1, "cm"),   # Space between columns
    legend.margin = margin(3, 3, 3, 3))

FLCpie  

