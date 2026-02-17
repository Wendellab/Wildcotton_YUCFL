#######################################################################
#######################################################################
library(ape)
library(ggtree)
library(treeio)
library(tidyverse)

YUCCultivar <- read.nexus("YUCCultivar_n277nj.tre")

YUCCultivar2 <- as_tibble(YUCCultivar) %>%
  mutate(label = gsub("SiPr_W9", "AD1_YUC_SiPr_W9", label),
         label = gsub("AD1_*", "", label),
         Population = str_extract(label, "[^_]*_[^_]*"),
         Population = case_when(
           grepl("^GD_", Population) ~ "GD",
           grepl("^Cultivar_", Population) ~ "Cultivar",
           grepl("^LR1_", Population) ~ "LR1",
           grepl("^LR2_", Population) ~ "LR2",
           TRUE ~ Population),
         Group = str_extract(label, "[^_]*")) %>%
  as.treedata()
  
  


YUCnjtree <- ggtree(YUCCultivar2, layout = "equal_angle") +
  geom_tippoint(aes(color=Group, shape=Group), size=4) +
  theme(legend.position =  c(0.7, 0.2),
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank()) # get rid of minor grid 

YUCnjtree

#######################################################################
### PCA ###############################################################
#######################################################################

pca <- read_table("YUCCultivar_n227.eigenvec", col_names = FALSE)
eigenval <- scan("YUCCultivar_n227.eigenval")

pca <- pca[,-1]
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pca2 <- pca %>% 
  mutate(ind = gsub("SiPr_W9", "AD1_YUC_SiPr_W9", ind),
         Population = gsub("AD1_*", "", ind),
         Population = str_extract(Population, "[^_]*_[^_]*"),
         Population = case_when(
           grepl("^GD_", Population) ~ "GD",
           grepl("^Cultivar_", Population) ~ "Cultivar",
           grepl("^LR1_", Population) ~ "LR1",
           grepl("^LR2_", Population) ~ "LR2",
           TRUE ~ Population),
         Population2 = str_extract(Population, "[^_]*"))


unique(pca2$Population2)

pca_plot <- ggplot(pca2,aes(PC1, PC2, color = Population2, shape = Population2)) + 
  geom_point(size=4) + 
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +  
  theme(axis.title=element_text(size=10),
        legend.position = "none",
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_))


structureplot <- YUCnjtree/pca_plot
structureplot
