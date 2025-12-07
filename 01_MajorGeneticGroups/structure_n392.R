setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig1_nj_structure_he_plastome")
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

YUC_distree <- read.newick("n392_nj_reroot.tre")

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
  mutate(Group = stringr::str_extract(Population, "[^_]*")) %>%
  mutate_at('Population', ~replace_na(.,"")) %>%
  as.treedata()
  

#YUC_distree3 = full_join(YUC_distree2, het, by = c('label' = 'INDV')) %>%
#  mutate(het = percent(het, accuracy = 0.1))   %>%
#  mutate(label2_het = paste(label2, het, sep = " ----- ")) %>%
#  as.treedata()

YUCnjtree <- ggtree(YUC_distree2, layout = "fan", show.legend = F) +
  geom_tippoint(aes(color=Population, shape=Group), size=4) +
  geom_rootedge(rootedge = 0.005, linewidth = 0.9) + 
  geom_tiplab(aes(label= label2), align = T, offset=0.005, size = 2, show.legend = F) +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values=safe_shape_palette) +
  xlim_tree(0.001) +
  theme(legend.position =  "none",
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank()) 


YUCnjtree

#######################################################################
#######################################################################

load("LEAstruture_n65.RData")

colplot <- plot_list[[28]]$data
colnames(colplot)[4] = 'pop_bar'
cbbPalette <- colorRampPalette( brewer.pal(n = 12, name = "Paired") )(28)

YUCnjtree2 <- YUCnjtree +  new_scale_fill() + 
  geom_fruit(data=colplot, geom=geom_col,
             mapping= aes(y=individual, x=q, fill=pop_bar),
             axis.params=list( axis = "x", text  = "M",  text.size  = 3, limits = c(0, 1.1)),
             offset = 0.3, pwidth = 0.5) +
  scale_fill_manual(values = cbbPalette) +
  theme(legend.position =  "none",
        panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.title.y = element_text(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        rect = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent",  colour = NA_character_),
        panel.background = element_rect(fill = "transparent", colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank())


#colplot2 <- plot_list[[29]]$data
#colnames(colplot2)[4] = 'pop_bar'

#YUCnjtree3 <- YUCnjtree2 +  new_scale_fill() + 
#  geom_fruit(data=colplot2, geom=geom_col,
#             mapping= aes(y=individual, x=q, fill=pop_bar),
#             axis.params=list( axis = "x", text  = "M",  text.size  = 3, limits = c(0, 1.2)),
#             offset = 0.02, pwidth = 0.8) +
#  scale_fill_manual(values = cbbPalette) +
#  theme(legend.position =  "none",
#        panel.spacing.x = unit(0, "lines"),
#        axis.line = element_blank(),
#        axis.title.y = element_text(),
#        axis.text.x = element_blank(),
#        strip.background = element_blank(),
#        panel.background = element_blank(),
#        axis.title = element_blank(),
#        panel.grid = element_blank())

YUCnjtree2$labels$colour <- "Group"
YUCnjtree2$labels$shape <- "Group"
YUCnjtree2$labels$fill <- "Population"

YUCnjtree2
#dev.off()


######################################################################
######################################################################

sampletable <- read.csv("samplegroupinfor.csv") %>%
  mutate(ID = gsub("_plastome.*", "", ID))

rownames(sampletable) <- sampletable$ID
mycha.melt <- melt(as.matrix(sampletable))
mycha.melt.group <- mycha.melt[mycha.melt$Var2 != 'ID',]  
colnames(mycha.melt.group)[3] <- "PlastomeType"
unique(mycha.melt.group$PlastomeType)

mycha.melt.group$PlastomeType <- gsub("AD2", "G. barbadense", mycha.melt.group$PlastomeType)
mycha.melt.group$PlastomeType <- gsub("Domesticated", "G. hirsutum Domesticated", mycha.melt.group$PlastomeType)
mycha.melt.group$PlastomeType <- gsub("Carribean ", "G. hirsutum Caribbean", mycha.melt.group$PlastomeType)
mycha.melt.group$PlastomeType <- gsub("Outgroup", "G. mustelinum", mycha.melt.group$PlastomeType)
mycha.melt.group$PlastomeType <- gsub("Yucatan", "G. hirsutum Yucatan", mycha.melt.group$PlastomeType)


###########


YUCnjtree3 <- YUCnjtree2 +  new_scale_fill() + 
  geom_fruit(data=mycha.melt.group, geom=geom_tile,
             mapping= aes(y=Var1, x=Var2, fill=PlastomeType),
             alpha = 0.8, 
             offset = 0.01, width = 0.01) +
  scale_fill_manual(values= c("G. barbadense"= "#31688EFF",
                              "G. hirsutum Yucatan" = "#CC79A7",
                              "G. hirsutum Domesticated" = "#35B779FF",
                              "G. hirsutum Caribbean" = "#D55E00",
                              "G. mustelinum" = "#FDE725FF")) +
  theme(legend.position = c(0.97,0.5),
        panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.title.y = element_text(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())



#######################################################################
#######################################################################




pdf("FigS1_YUCFL_NJ_Structure_Plastome_n392_FINAL.pdf", width = 23, height = 22)
YUCnjtree3
dev.off()











#######################################################################
### PCA ###############################################################
#######################################################################

cbbPalette2 <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", 
                 "#0072B2", "#D55E00",  "#999999", "#000000")

pca <- read_table("YUCFLAD2AD4_n392.eigenvec", col_names = FALSE)
eigenval <- scan("YUCFLAD2AD4_n392.eigenval")

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
  mutate(Population = gsub("AD1_*", "", Population)) %>%
  mutate(Population = stringr::str_extract(Population, "[^_]*")) %>%
  mutate(Population = gsub("GD_.*", "GD", Population)) %>%
  mutate(Population = gsub("Cultivar_.*", "Cultivar", Population)) %>%
  mutate(Population = gsub("LR1_.*", "LR1", Population)) %>%
  mutate(Population = gsub("LR2_.*", "LR2", Population)) %>%
  mutate_at('Population', ~replace_na(.,"")) 


unique(pca2$Population)

pca_plot <- ggplot(pca2,aes(PC1, PC2, fill = Population, color = Population, shape = Population)) + 
  geom_point(size=4) + 
  #  scale_fill_manual(values= safe_colorblind_palette, guide = "none") +
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values = safe_shape_palette) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +  
  theme(legend.position =  c(0.1,0.8),
        legend.background = element_blank(),
        legend.box = "horizontal", 
        #legend.box.background = element_rect(linetype="dotted", colour = "black"),
        #legend.title=element_blank(),
        legend.key = element_rect(fill = "white", colour = "grey30", linetype="dotted")
  ) 

pca_plot




pdf("FigS0_YUCFL_n392.pdf", width = 10, height = 10)
pca_plot
dev.off()












#######################################################################
#######################################################################
cbbPalette3 <- colorRampPalette( brewer.pal(n = 12, name = "Paired") )(12)


load("LEAstruture_FLn166.RData")
colplot3 <- plot_list[[12]]$data
colnames(colplot3)[4] = 'pop_bar'

load("LEAstruture_YUCn158.RData")
colplot4 <- plot_list[[10]]$data
colnames(colplot4)[4] = 'pop_bar'



YUCnjtree5 <- YUCnjtree4 +  new_scale_fill() + 
  geom_fruit(data=colplot3, geom=geom_col,
             mapping= aes(y=individual, x=q, fill=pop_bar),
             axis.params=list( axis = "x", text  = "M",  text.size  = 3, limits = c(0, 1.2)),
             offset = 0.2, pwidth = 0.8) +
  scale_fill_manual(values = cbbPalette3) +
  
  new_scale_fill() + 
  geom_fruit(data=colplot4, geom=geom_col,
             mapping= aes(y=individual, x=q, fill=pop_bar),
             axis.params=list( axis = "x", text  = "M",  text.size  = 3, limits = c(0, 1.2)),
             offset = -0.8, pwidth = 0.8) +
  scale_fill_manual(values = cbbPalette3) +
  
  theme(legend.position =  "none",
        panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.title.y = element_text(),
        axis.text.x = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank())

YUCnjtree5




pdf("FigS2_YUCFL_n392.pdf", width = 13, height = 30)
YUCnjtree5
dev.off()










##############################################################
#############relatedness######################################
##############################################################


# https://stackoverflow.com/questions/9617348/reshape-three-column-data-frame-to-matrix-long-to-wide-format
library(dplyr)
library("stringr") 

relatFL <- read.table("FLn166.genome", header = T) %>%
  select(IID1, IID2, PI_HAT) %>%
  mutate(IID1 = gsub("AD1_", "", IID1)) %>%
  mutate(IID2 = gsub("AD1_", "", IID2)) %>%
  # mutate(PI_HAT = if_else(PI_HAT < 0.3, NA, PI_HAT)) %>%
  
  # filter( grepl("PR325",IID1)) %>%
  # filter( grepl("PR325",IID2)) %>%
  # filter( !grepl("GD2",IID2)) %>%
  mutate(IID2 = fct_rev(fct_infreq(IID2)),
         IID1 = fct_infreq(IID1))

mean(relatFL$PI_HAT, na.rm = T)
subset(relatFL,PI_HAT != 0)


relatFLplot <- ggplot(data=relatFL_lower, aes(IID2, IID1, fill = PI_HAT))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white",
                       high = "red",
                       mid = "pink",
                       midpoint = 0.5, 
                       limit = c(0,1),
                       space = "Lab") +
  geom_text(data=subset(relatFL,PI_HAT > 0.3), aes(label = round(PI_HAT,2)), size = 2.5) +
  #  theme_minimal()+ 
  theme(axis.text.x =element_text(size=13,  angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  coord_fixed()

relatFLplot


##################


relatYUC <- read.table("YUCn158.genome", header = T) %>%
  select(IID1, IID2, PI_HAT) %>%
  mutate(IID1 = gsub("AD1_", "", IID1)) %>%
  mutate(IID2 = gsub("AD1_", "", IID2)) %>%
  # mutate(PI_HAT = if_else(PI_HAT < 0.3, NA, PI_HAT)) %>%
  
  # filter( grepl("PR325",IID1)) %>%
  # filter( grepl("PR325",IID2)) %>%
  # filter( !grepl("GD2",IID2)) %>%
  mutate(IID2 = fct_rev(fct_infreq(IID2)),
         IID1 = fct_infreq(IID1))

mean(relatYUC$PI_HAT, na.rm = T)

relatYUCplot <- ggplot(data=subset(relatYUC,PI_HAT != 0), aes(IID1, IID2, fill = PI_HAT))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white",
                       high = "red",
                       mid = "pink",
                       midpoint = 0.5, 
                       limit = c(0,1),
                       space = "Lab") +
  geom_text(data=subset(relatYUC,PI_HAT > 0.3), aes(label = round(PI_HAT,2)), size = 2.5) +
  #  theme_minimal()+ 
  theme(axis.text.x =element_text(size=13,  angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size=13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  coord_fixed()

relatYUCplot

pdf("FigS3_YUCFL_n392.pdf", width = 30, height = 30)
relatFLplot
relatYUCplot
dev.off()

