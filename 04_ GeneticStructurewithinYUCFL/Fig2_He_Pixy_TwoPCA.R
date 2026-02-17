setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig5_PopGeneStructure_YUCFL/Fig5_He_Pixy_Popstructure/")

library(dplyr)
library(ggplot2)
library(directlabels)
library(stringr)
library(geodist)
library(ape)
library(vegan)
library(smplot2)
library(tibble)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
#####################################################################
####################################################################

FLcbbPalette <- colorRampPalette(brewer.pal(n = 12, name = "Paired"))(15)
FLcategories <- c("FL_CPH", "FL_CPT", "FL_FMY", "FL_NP", "FL_OPB2", 
                  "FL_OPB4", "FL_OPB5", "FL_PK", "FL_RBD", "FL_RBT", 
                  "FL_RNRB", "FL_SR", "FL_TC", "FL_VKPA", "FL_MK")
FLcolormapping <- setNames(FLcbbPalette, FLcategories)

YUCcbbPalette <- colorRampPalette(brewer.pal(n = 12, name = "Paired"))(8)
YUCcategories <- c("YUC_CeCo", "YUC_CeDo", "YUC_CeDr", "YUC_CePr", 
                   "YUC_RiCa", "YUC_RiCh", "YUC_SiPa", "YUC_SiPr")
YUCcolormapping <- setNames(YUCcbbPalette, YUCcategories)

safe_colorblind_palette <- c("LR1" = "#F4A582",
                             "LR2" = "#D6604D",
                             "Cultivar" = "#B2182B",
                             "PR"= "black",
                             "GD" = "#4D9221",
                             "FL" = "#3B8ABE",
                             "YUC" = "#8C66AF",
                             FLcolormapping,
                             YUCcolormapping)



safe_shape_palette <- c("FL"= 19,
                        "LR1" = 17,
                        "Cultivar" = 17,
                        "LR2" = 17,
                        "GD" = 3,
                        "PR" = 13,
                        "YUC" = 15,
                        'AD2'= 1,
                        'AD4'= 2)

#####################################################################
####################################################################
####################################################################
gpsmk <- read.csv("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig0_map/CottonSamplingMap5_ut8.csv", encoding="UTF-8")
gpsmk3 <- gpsmk %>%
  filter(Source == "Field") %>%
  filter(Region %in% c("Yucatan", "Florida"))  %>%
  filter(Group != "Landrace") %>%
  select(Region, ind, Pop_Site, lat, long )%>%
  mutate(ind = gsub("í","i", ind)) %>%
  mutate(Pop_Site = gsub("AD1_MK","AD1_FL_MK", Pop_Site)) %>%
  mutate(Individual = c(Pop_Site[1:166], ind[167:324])) %>%
  mutate(Site = str_extract(Individual, "[^_]*_[^_]*_[^_]*"))%>%
  mutate(Site = gsub("AD1_","", Site)) %>%
  select(Region, Individual, Site, lat, long )%>%
  add_count(Site, name = 'id_occurence') 

gpsmk4 <- gpsmk3 %>%
  select(!Individual) %>%
  mutate(Site2 = str_extract(Site, "_[^_]*"))%>%
  mutate(Site2 = gsub("_","", Site2)) %>%
  mutate(group = Site) %>%
  distinct(Site, .keep_all = TRUE)


#####################################################################
####################################################################
####################################################################


FLgeoDxy <- gpsmk4 %>%
  filter(Region == "Florida") %>%
  select("Site2", "lat", "long")  %>%
  column_to_rownames(var = "Site2") 

FLgeodis <- as.matrix(geodist(FLgeoDxy[,c("lat","long")], measure = "haversine"))
rownames(FLgeodis) <- row.names(FLgeoDxy)
colnames(FLgeodis) <- row.names(FLgeoDxy)

FLgeodis.melt <- melt(replace(FLgeodis, lower.tri(FLgeodis, TRUE), NA), na.rm = TRUE)%>%
  #mutate(group = paste(Var1, Var2, sep = "_")) 
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         group = ifelse(Var1 < Var2, paste(Var1, Var2, sep = "_"), paste(Var2, Var1, sep = "_"))) %>%
  rename("geodis" = "value") %>%
  select("group", "geodis")

FLdxy <- read.table("AD1350_25Pop.dxy4.txt", sep = "") %>%
  filter(startsWith(V1, "FL_") & startsWith(V2, "FL_")) %>%
  mutate( V1 = gsub("^FL_", "", V1),V2 = gsub("^FL_", "", V2))%>%
  mutate(V1 = as.character(V1),
         V2 = as.character(V2),
         group = ifelse(V1 < V2, paste(V1, V2, sep = "_"), paste(V2, V1, sep = "_")))%>% 
  rename("genedis" = "V3") %>%
  select("group", "genedis")

FLdxygeo <- merge(FLdxy,FLgeodis.melt, by = "group", all = T)


FLdxygeo.plot  <- ggplot(data = FLdxygeo, mapping = aes(x = geodis/1000, y = genedis)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
  # geom_text_repel(aes(label=group), max.overlaps=Inf, size = 2, show.legend = F) +
  sm_statCorr(color = "#0f993d", corr_method = "pearson") +
  xlab("FL Distance (km)") + ylab("FL Dxy") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x =element_text(size=10),
        axis.title.y =element_text(size=10),
        panel.border = element_rect(colour = "black", fill=NA))



FLdxygeo.plot
####################################################################
####################################################################
####################################################################

YUCgeoDxy <- gpsmk4 %>%
  filter(Region == "Yucatan") %>%
  select("Site2", "lat", "long")  %>%
  column_to_rownames(var = "Site2") 

YUCgeodis <- as.matrix(geodist(YUCgeoDxy[,c("lat","long")], measure = "haversine"))
rownames(YUCgeodis) <- row.names(YUCgeoDxy)
colnames(YUCgeodis) <- row.names(YUCgeoDxy)

YUCgeodis.melt <- melt(replace(YUCgeodis, lower.tri(YUCgeodis, TRUE), NA), na.rm = TRUE)%>%
  mutate(Var1 = as.character(Var1),
         Var2 = as.character(Var2),
         group = ifelse(Var1 < Var2, paste(Var1, Var2, sep = "_"), paste(Var2, Var1, sep = "_"))) %>%
  rename("geodis" = "value") %>%
  select("group", "geodis")

YUCdxy <- read.table("AD1350_25Pop.dxy4.txt") %>%
  filter(startsWith(V1, "YUC_") & startsWith(V2, "YUC_")) %>%
  mutate( V1 = gsub("^YUC_", "", V1),V2 = gsub("^YUC_", "", V2))%>%
  mutate(V1 = as.character(V1),
         V2 = as.character(V2),
         group = ifelse(V1 < V2, paste(V1, V2, sep = "_"), paste(V2, V1, sep = "_")))%>% 
  rename("genedis" = "V3") %>%
  select("group", "genedis")

YUCdxygeo <- merge(YUCdxy,YUCgeodis.melt, by = "group", all = T)


YUCdxygeo.plot  <- ggplot(data = YUCdxygeo, mapping = aes(x = geodis/1000, y = genedis)) +
  geom_point(shape = 21, fill = '#0f993d', color = 'white', size = 3) +
  # geom_text_repel(aes(label=group), max.overlaps=Inf, size = 2, show.legend = F) +
  sm_statCorr(color = "#0f993d", corr_method = "pearson")+
  xlab("YUC Distance (km)") + ylab("YUC Dxy") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x =element_text(size=10),
        axis.title.y =element_text(size=10),
        panel.border = element_rect(colour = "black", fill=NA))


########################################################################
########################################################################





library(gridExtra)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(grid)
library(ggplot2)
library(ape)
library(ggtree)
library(treeio)
library(rcartocolor)
library(RColorBrewer)
library(ggnewscale)
library(ggtreeExtra)
library(reshape2)
library(scales)

#######################################################################
### PCA ###############################################################
#######################################################################

pca <- read_table("FLn166.eigenvec", col_names = FALSE)
eigenval <- scan("FLn166.eigenval")

pca <- pca[,-1]
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

names(pca)[1] <- "ind"
pca$ind <- gsub("Pop","Site",pca$ind)

# set names
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pca2 <- pca %>% 
  mutate(Population = gsub("_peru", "", ind)) %>%
  #mutate(label2 = gsub("AD2_*", "Gb_", label2)) %>%
  #mutate(label2 = gsub("AD4_mus*", "Gm_outgroup", label2)) %>%
  mutate(Population = gsub("AD1_*", "", Population)) %>%
  mutate(Population = stringr::str_extract(Population, "[^_]*_[^_]*")) %>%
  mutate(Population2 = stringr::str_extract(Population, "[^_]*")) %>%
  mutate(Population3 = gsub("FL_", "", Population)) %>%
  mutate_at('Population', ~replace_na(.,""))  %>%
  mutate_at('Population2', ~replace_na(.,""))


unique(pca2$Population2)

pca_plotFL <- ggplot(pca2,aes(PC1, PC2, color = Population, shape = Population2)) + 
  geom_point(size=4) + 
  geom_text_repel( data=subset(pca2, Population == "FL_RNRB" | !duplicated(pca2$Population)),
                   aes(label=Population3, color = Population), max.overlaps=Inf, 
                   size = 3, show.legend = F) +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values = safe_shape_palette) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +  
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x =element_text(size=10),
        axis.title.y =element_text(size=10),
        legend.position = "none",
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_)) # necessary to avoid drawing panel outline
#panel.grid.major = element_blank(), # get rid of major grid
#panel.grid.minor = element_blank()) # get rid of minor grid 


#######################################################################
### PCA ###############################################################
#######################################################################

pca <- read_table("YUCn158.eigenvec", col_names = FALSE)
eigenval <- scan("YUCn158.eigenval")

pca <- pca[,-1]
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

names(pca)[1] <- "ind"
pca$ind <- gsub("Pop","Site",pca$ind)

# set names
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pca2 <- pca %>% 
  mutate(Population = gsub("_peru", "", ind)) %>%
  #mutate(label2 = gsub("AD2_*", "Gb_", label2)) %>%
  #mutate(label2 = gsub("AD4_mus*", "Gm_outgroup", label2)) %>%
  mutate(Population = gsub("AD1_*", "", Population)) %>%
  mutate(Population = stringr::str_extract(Population, "[^_]*_[^_]*")) %>%
  mutate(Population2 = stringr::str_extract(Population, "[^_]*")) %>%
  mutate(Population3 = gsub("YUC_", "", Population)) %>%
  mutate_at('Population', ~replace_na(.,""))  %>%
  mutate_at('Population2', ~replace_na(.,""))


unique(pca2$Population2)

pca_plotYUC <- ggplot(pca2,aes(PC1, PC2, color = Population, shape = Population2)) + 
  geom_point(size=4) + 
  geom_text_repel( data=subset(pca2, !duplicated(pca2$Population)),
                   aes(label=Population3, color = Population), max.overlaps=Inf, 
                   size = 3, show.legend = F) +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values = safe_shape_palette) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +  
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x =element_text(size=10),
        axis.title.y =element_text(size=10),
        legend.position = "none",
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_)) # necessary to avoid drawing panel outline



########################################################################
########################################################################
########################################################################
library(cowplot)

finalplot <- ggdraw() +
  draw_plot(pca_plotYUC, x=0, y=0.5, width=1/3, height=0.5) +
  draw_plot(YUCpie, x=1/3, y=0.5, width=1/3, height=0.5) +
  draw_plot(YUCdxygeo.plot, x=2/3, y=0.5, width=1/3, height=0.5) +
  draw_plot(pca_plotFL, x=0, y=0, width=1/3, height=0.5) +
  draw_plot(FLCpie, x=1/3, y=0, width=1/3, height=0.5) +
  draw_plot(FLdxygeo.plot, x=2/3, y=0, width=1/3, height=0.5) +
  draw_plot_label(label=c("a","b","c","d","e","f"), size=16, fontface="bold",
                  x=c(0,1/3,2/3,0,1/3,2/3), y=c(1,1,1,0.5,0.5,0.5))

pdf("../Fig5_Het_PopStructure_GeovsGenetic_FLYUC_v2.pdf", width = 16, height = 10)
finalplot
dev.off()
