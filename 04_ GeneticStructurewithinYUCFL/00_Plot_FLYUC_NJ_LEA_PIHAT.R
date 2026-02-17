setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig5_PopGeneStructure_YUCFL/FigS678_nj_structure_relatedness/")
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
library(ggbreak)

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


###########################YUC samples below###########################
#######################################################################

YUC_distree <- read.newick("YUCn188_nj_rt.tre")


YUC_distree2 <- as_tibble(YUC_distree) %>% 
  mutate(label2 = gsub("_peru", "", label)) %>%
  #mutate(label2 = gsub("AD2_*", "Gb_", label2)) %>%
  #mutate(label2 = gsub("AD4_mus*", "Gm_outgroup", label2)) %>%
  mutate(Population = gsub("AD1_*", "", label2)) %>%
  mutate(Population = stringr::str_extract(Population, "[^_]*_[^_]*")) %>%
  mutate(Population = gsub("GD_.*", "GD", Population)) %>%
  mutate(Population = gsub("Cultivar_.*", "Cultivar", Population)) %>%
  mutate(Population = gsub("LR1_.*", "LR1", Population)) %>%
  mutate(Population = gsub("LR2_.*", "LR2", Population)) %>%
  mutate(label2 = gsub("AD1_*", "", label2)) %>%
  mutate(label = gsub("AD1_*", "", label)) %>%
  mutate(Group = stringr::str_extract(Population, "[^_]*")) %>%
  mutate_at('Population', ~replace_na(.,"")) %>%
  as.treedata()


YUCnjtree <- ggtree(YUC_distree2,  show.legend = F) +
  geom_tippoint(aes(color=Population, shape = Group), size=2) +
  geom_rootedge(rootedge = 0.005, linewidth = 0.9) + 
  geom_tiplab(aes(label= label2), align = T, offset=0.005, size = 2, show.legend = F) +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values=safe_shape_palette) +
  xlim_tree(0.5) +
  theme(legend.position =  "none") + ylim(-10, 190) 

YUCnjtree
#YUCnjtree <- YUCnjtree0 %>% rotate(196) %>% rotate(277)  


#p <- ggtree(YUC_distree2) + geom_text(aes(label=node)) 
#p  %>% rotate(282) %>% rotate(283)

#######################################################################
#######################################################################

load("YUCn188_LEAstruture_n65.RData")

colplot <- plot_list[[11]]$data %>%
  mutate(individual = gsub("AD1_*", "", individual)) 
colnames(colplot)[4] = 'pop_bar'
cbbPalette <- colorRampPalette(brewer.pal(n = 12, name = "Paired"))(12)
colplot <- colplot %>%
  mutate(pop_bar = factor(pop_bar, levels = paste0("P", 1:12)))

YUCnjtree2 <- YUCnjtree +  new_scale_fill() + 
  geom_fruit(data=colplot, geom=geom_col,
             mapping= aes(y=individual, x=q, fill=pop_bar),
             axis.params=list( axis = "x", text  = "M",  text.size  = 3, limits = c(0, 1.2)),
             offset = 0.33, pwidth = 0.5) +
  scale_fill_manual(values = YUC_component_color_mapping) +
  theme(legend.position =  "none",
        panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.title.y = element_text(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

YUCnjtree2


#######################################################################
#######################################################################
YUCfactor_data <- as.factor(gsub("\n", "", unlist(YUC_distree2@phylo$tip.label)))

relatYUC <- read.table("YUCn158.genome", header = T) %>%
  select(IID1, IID2, PI_HAT) %>%
  mutate(IID1 = gsub("AD1_", "", IID1)) %>%
  mutate(IID2 = gsub("AD1_", "", IID2)) %>%
  mutate(across(PI_HAT, as.numeric)) %>%
  filter(PI_HAT > 0.2) %>%
  mutate(IID2 = factor(IID2, levels = intersect(YUCfactor_data, IID2)))
  #mutate(PI_HAT = if_else(PI_HAT < 0.3, NA, PI_HAT)) 

relatYUC %>% 
  filter(grepl("RiCh|RiCa", IID1), grepl("RiCh|RiCa", IID2), as.numeric(PI_HAT) > 0) %>% 
  summarise(mean_PI_HAT = mean(as.numeric(PI_HAT), na.rm = TRUE))


relatYUC %>%
  filter(!(grepl("RiCh|RiCa", IID1) & grepl("RiCh|RiCa", IID2)),
         as.numeric(PI_HAT) > 0) %>%
  summarise(mean_PI_HAT = mean(as.numeric(PI_HAT), na.rm = TRUE))


YUCnjtree3 <- YUCnjtree2  +  new_scale_fill() + 
  geom_fruit(data=relatYUC, geom=geom_tile,
             mapping= aes( y=IID1, x= IID2, fill = PI_HAT),
             axis.params=list( axis = "x", text  = "M",  text.size  = 2, text.angle = -90, hjust = 0.1),
             offset = 0.01, pwidth = 3.5) +
  scale_fill_gradientn( colors = c("white", "lightblue", "yellow", "orange", "red"),  # More distinct colors
                        values = scales::rescale(c(0.2, 0.3, 0.6, 1)),                 # Breaks at 0.3 and 0.6
                        limits = c(0.2, 1),
                        name = "Scale") +
  theme_void() +
  theme(legend.position = "right",
        legend.justification = c(0, 0.5), # anchor right-center
        axis.title.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank()) 

YUCnjtree3

pdf("FigS7_YUCn158_nj_structure_PIhat.pdf", width = 18, height = 13)
YUCnjtree3
dev.off()

############################################################################
############################################################################





##########################FL samples below##################################
############################################################################
FL_distree <- read.newick("FLn196_nj_rt.tre")


FL_distree2 <- as_tibble(FL_distree) %>% 
  mutate(label2 = gsub("_peru", "", label)) %>%
  #mutate(label2 = gsub("AD2_*", "Gb_", label2)) %>%
  #mutate(label2 = gsub("AD4_mus*", "Gm_outgroup", label2)) %>%
  mutate(Population = gsub("AD1_*", "", label2)) %>%
  mutate(Population = stringr::str_extract(Population, "[^_]*_[^_]*")) %>%
  mutate(Population = gsub("GD_.*", "GD", Population)) %>%
  mutate(Population = gsub("Cultivar_.*", "Cultivar", Population)) %>%
  mutate(Population = gsub("LR1_.*", "LR1", Population)) %>%
  mutate(Population = gsub("LR2_.*", "LR2", Population)) %>%
  mutate(label2 = gsub("AD1_*", "", label2)) %>%
  mutate(label = gsub("AD1_*", "", label)) %>%
  mutate(Group = stringr::str_extract(Population, "[^_]*")) %>%
  mutate_at('Population', ~replace_na(.,"")) %>%
  as.treedata()


FLnjtree <- ggtree(FL_distree2,  show.legend = F) +
  geom_tippoint(aes(color=Population, shape = Group), size=2, guide="none") +
  geom_rootedge(rootedge = 0.005, linewidth = 0.9) + 
  geom_tiplab(aes(label= label2), align = T, offset=0.005, size = 2, show.legend = F) +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values=safe_shape_palette) +
  xlim_tree(0.5) +
  theme(legend.position =  "none")  + 
  #scale_x_break(c(0.002, 0.087)) + 
  ylim(-3, 200) 

#######################################################################
#######################################################################

load("FLn196_LEAstruture.RData")

colplot <- plot_list[[15]]$data %>%
  mutate(individual = gsub("AD1_*", "", individual)) 
colnames(colplot)[4] = 'pop_bar'
cbbPalette <- colorRampPalette(brewer.pal(n = 12, name = "Paired"))(15)
# Convert pop_bar to numeric-sorted factor
colplot <- colplot %>%
  mutate(pop_bar = factor(pop_bar, levels = paste0("P", 1:15)))


levels(as.factor(colplot$pop_bar))

FLnjtree2 <- FLnjtree +  new_scale_fill() + 
  geom_fruit(data=colplot, geom=geom_col,
             mapping= aes(y=individual, x=q, fill=pop_bar),
             axis.params=list( axis = "x",  text.size  = 3, limits = c(0, 1.2)),
             offset = 0.3, pwidth = 0.5) +
  scale_fill_manual(values = FL_component_color_mapping) +
  theme(legend.position =  "none",
        panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.title.y = element_text(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

FLnjtree2


#######################################################################
#######################################################################
FLfactor_data <- as.factor(gsub("\n", "", unlist(FL_distree2@phylo$tip.label)))

relatFL <- read.table("FLn166.genome", header = T) %>%
  select(IID1, IID2, PI_HAT) %>%
  mutate(IID1 = gsub("AD1_", "", IID1)) %>%
  mutate(IID2 = gsub("AD1_", "", IID2)) %>%
  mutate(across(PI_HAT, as.numeric)) %>%
  filter(PI_HAT > 0.2) %>%
  mutate(IID2 = factor(IID2, levels = intersect(FLfactor_data, IID2)))

relatFL %>%
  filter(
    grepl("MK|OPB2|SR|NP|TC", IID1),
    grepl("MK|OPB2|SR|NP|TC", IID2),
    as.numeric(PI_HAT) > 0
  ) %>%
  summarise(mean_PI_HAT = mean(as.numeric(PI_HAT), na.rm = TRUE))

relatFL %>%
  filter(
    #!grepl("MK|OPB2|SR|NP|TC", IID1) & !grepl("MK|OPB2|SR|NP|TC", IID2),
    as.numeric(PI_HAT) > 0
  ) %>%
  summarise(mean_PI_HAT = mean(as.numeric(PI_HAT), na.rm = TRUE))


FLnjtree3 <- FLnjtree2  +  new_scale_fill() + 
  geom_fruit(data=relatFL, geom=geom_tile,
             mapping= aes( y=IID1, x= IID2, fill = PI_HAT),
             axis.params=list( axis = "x", text.size  = 2, text.angle = -90, hjust = 0.1),
             offset = 0.01, pwidth = 2.5) +
  scale_fill_gradientn( colors = c("white", "lightblue", "yellow", "orange", "red"),  # More distinct colors
                        values = scales::rescale(c(0.2, 0.3, 0.6, 1)),                 # Breaks at 0.3 and 0.6
                        limits = c(0, 1),
                        name = "Scale") +
  theme_void() +
  theme(legend.position = c(0.06,0.65),
        axis.title.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank()) 


FLnjtree3


#######################################################################
#######################################################################

pdf("FigS8_FLn166_nj_structure_PIhat.pdf", width = 20, height = 15)
FLnjtree3 + guides(#fill=guide_legend(ncol=2),
                   color = guide_legend(ncol=2),
                   shape=guide_legend(ncol=2))
dev.off()

