setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig1_nj_structure_he_plastome/Fig1PCA/")
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

#####################################################################
#####################################################################
#FLcolors <- colorRampPalette(c("#014636", "#01665E", "#35978F", "#80CDC1", "#C7EAE5"))(15)

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


#######################################################################
#######################################################################

YUC_distree <- read.nexus("nj_n380.tre")

YUC_distree2 <- as_tibble(YUC_distree) %>%
  mutate(
    label2 = gsub("_peru|AD1_*", "", label),
    Population = str_extract(label2, "[^_]*_[^_]*"),
    Population = case_when(
      grepl("^GD_", Population) ~ "GD",
      grepl("^Cultivar_", Population) ~ "Cultivar",
      grepl("^LR1_", Population) ~ "LR1",
      grepl("^LR2_", Population) ~ "LR2",
      TRUE ~ Population
    ),
    Group  = sapply(strsplit(Population, "_"), `[`, 1),
    Group2 = sapply(strsplit(Population, "_"), function(x) if(length(x) > 1) x[2] else x[1])
  ) %>%
  mutate(across(c(Population, Group, Group2), ~replace_na(.,"")))

labeltree <- YUC_distree2 %>% 
  filter(!is.na(Group2)) %>%
  filter(!grepl("Ph", Group2)) %>%     
  distinct(Group2, .keep_all = TRUE) %>%
  select(node) %>%
  pull() 


YUC_distree2 <- as.treedata(YUC_distree2)

YUC_distree2@data$Group <- factor(as.factor(YUC_distree2@data$Group), 
                     levels = c("AD4","AD2","Cultivar", "LR1","LR2","PR", "GD", "FL", "YUC")) 


YUCnjtree <- ggtree(YUC_distree2, layout = "equal_angle") +
  geom_tippoint(aes(color=Population, shape=Group), size=4) +
  geom_rootedge(rootedge = 0.005, linewidth = 0.9) + 
  geom_tiplab(aes(subset=(node %in% labeltree), label= Group2, color = Population),size = 4, hjust = -1) +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values=safe_shape_palette) +
  theme(legend.position =  "none",
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank()) # get rid of minor grid 

YUCnjtree

dev.off()


#######################################################################
### PCA ###############################################################
#######################################################################

pca <- read_table("AD1380.eigenvec", col_names = FALSE)
eigenval <- scan("AD1380.eigenval")

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
  mutate(Population = gsub("GD_.*", "GD", Population)) %>%
  mutate(Population = gsub("Cultivar_.*", "Cultivar", Population)) %>%
  mutate(Population = gsub("LR1_.*", "LR1", Population)) %>%
  mutate(Population = gsub("LR2_.*", "LR2", Population)) %>%
  mutate(Population2 = stringr::str_extract(Population, "[^_]*")) %>%
  mutate_at('Population', ~replace_na(.,""))  %>%
  mutate_at('Population2', ~replace_na(.,""))


unique(pca2$Population2)

pca_plot <- ggplot(pca2,aes(PC1, PC2, color = Population, shape = Population2)) + 
  geom_point(size=4) + 
  #  scale_fill_manual(values= safe_colorblind_palette, guide = "none") +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values = safe_shape_palette) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +  
  theme(axis.title=element_text(size=14),
        legend.position = "none",
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_)) # necessary to avoid drawing panel outline
        #panel.grid.major = element_blank(), # get rid of major grid
        #panel.grid.minor = element_blank()) # get rid of minor grid 

pca_plot + guides(colour = "none") 

#######################################################################
#######################################################################


pca <- read_table("YUCFLAD2AD4_n392.eigenvec", col_names = FALSE)
eigenval <- scan("YUCFLAD2AD4_n392.eigenval")

pca <- pca[,-1]
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

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
  mutate(Population = gsub("GD_.*", "GD", Population)) %>%
  mutate(Population = gsub("Cultivar_.*", "Cultivar", Population)) %>%
  mutate(Population = gsub("LR1_.*", "LR1", Population)) %>%
  mutate(Population = gsub("LR2_.*", "LR2", Population)) %>%
  mutate(Population2 = stringr::str_extract(Population, "[^_]*")) %>%
  mutate_at('Population', ~replace_na(.,""))  %>%
  mutate_at('Population2', ~replace_na(.,""))


unique(pca2$Population2)

pca_plot2 <- ggplot(pca2,aes(PC1, PC2, color = Population, shape = Population2)) + 
  geom_point(size=4) + 
  scale_color_manual(values= safe_colorblind_palette, guide = "none") +
  scale_shape_manual(values = safe_shape_palette) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +  
  theme(legend.position = "none",
    legend.background = element_blank(),
    rect = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent",  colour = NA_character_),
    panel.background = element_rect(fill = "transparent", colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    #axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title=element_blank()) 


library(cowplot)

PCAinsert <- ggdraw() +
  draw_plot(pca_plot, 0, 0, 1, 1) +  # Main plot takes full area
  draw_plot(pca_plot2, 0.08, 0.68, 0.3, 0.3)  # Inset plot in top-right corner

PCAinsert

finalplot2 <- ggdraw() +
  draw_plot(mapfinalplot, x = 0, y = 0.5, width = 1, height = 0.5) +
  
  draw_plot(PCAinsert, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(YUCnjtree, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 20,  fontface = "bold",
                  x = c(0,0.58,0,0.5), y = c(1,1,0.5,0.5))


pdf("Fig1_PCA_NJtree_n380_v3.pdf", width = 16, height = 16)
finalplot2
dev.off()

#######################################################################
#######################################################################


