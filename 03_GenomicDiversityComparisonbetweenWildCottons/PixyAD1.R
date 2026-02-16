setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig3_Pixy_He_FIS_ROH/Fig3_Pixy_AD1n380/")
library(stringr)
library(dplyr)
library(grid)   # for the textGrob() function
library(ggpubr)
library(ggplot2)
library(scales)
library(reshape2)
library(pheatmap)
library("cowplot")
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(forcats)
#######################################################################


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
                             "Wild" = "#E69FAD",
                             "PR"= "black",
                             "GD" = "#4D9221",
                             "FL" = "#3B8ABE",
                             "YUC-E" = "#8C66AF",
                             "YUC-W" = "#8C66AF",
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
### Pi ###############################################################
#######################################################################

inp.pi <- read.table("AD1380_8groups.pi3.txt",sep = "",header=F) %>% 
  rename("pop" = V1,
         "chromosome" = V2,
         "avg_pi" = V3,
         "pi_sd" = V4)  %>%
  mutate(pop = recode(pop,"YUCE" = "YUC-E", "YUCW" = "YUC-W"),
         group = case_when(
           pop %in% c("Cultivar", "LR1", "LR2", "Wild") ~ "Global",
           TRUE ~ "Regional"))


inp.pi$pop <- gsub("AD1_", "", inp.pi$pop)
group_ordered <- with(inp.pi,  reorder(pop,  avg_pi,  mean))
inp.pi$pop <- factor(inp.pi$pop,levels = levels(group_ordered))

inp.pi$pop <- factor(inp.pi$pop,  levels = c("Cultivar", "LR1", "LR2","Wild", "GD", "PR","FL", "YUC-E", "YUC-W"))

inp.pi.avg.mean <- inp.pi %>%
  group_by(pop, group) %>% 
  dplyr::summarize(Avg = mean(avg_pi, na.rm=TRUE)) %>%
  mutate(label = gsub("e", " %*% 10^", sprintf("%.2e", Avg)))%>%
  as.data.frame()

plot.inp.pi.avg <- ggplot(inp.pi, aes(x = pop, y = avg_pi, fill = pop, color = pop)) +
  geom_violin(alpha = 0.3) +
  scale_fill_manual(values = safe_colorblind_palette) +
  #scale_shape_manual(values = safe_shape_palette) +
  facet_grid(. ~ group, scales = "free_x", space = "fixed")+
  geom_boxplot(width = 0.1, color = "white", alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.3, color = "black", alpha = 0.5, show.legend = FALSE) +  # dots
  ylab("Avg Pi / Chromosome") + 
  theme_minimal() +
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    legend.position="none",
    axis.text.x = element_text(size = 10, color = "black")
    # axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1, size = 4)
  )

plot.inp.pi.avg2 <- plot.inp.pi.avg +
  geom_text(data = inp.pi.avg.mean, aes(x = pop, y = 0, 
                                        label = gsub("e([+-]?)([0-9]+)", " %*% 10^\\1\\2", sprintf("%.2e", Avg))), 
            size =2.7,parse = TRUE,  color = "blue")


plot.inp.pi.avg2


##############################################
##############################################
# Read and process
inp.dxy <- read.table("AD1380_8groups.dxy3.txt", sep = "", header = FALSE) %>% 
  rename(pop1 = V1, pop2 = V2, chromosome = V3, avg_dxy = V4, dxy_sd = V5) %>%
  mutate(pop1 = recode(pop1, "YUCE" = "YUC-E", "YUCW" = "YUC-W"),
         pop2 = recode(pop2, "YUCE" = "YUC-E", "YUCW" = "YUC-W"))

pop_levels <- c("Cultivar", "LR1", "LR2", "GD", "PR","FL", "YUC-E", "YUC-W")

# Get averages
inp.avg <- inp.dxy %>%
  group_by(pop1, pop2) %>%
  summarise(Mean = mean(avg_dxy) * 1000, .groups = "drop") %>%
  mutate(pop1 = factor(pop1, pop_levels),
         pop2 = factor(pop2, pop_levels))

# Create full matrix and fill symmetric
full_grid <- expand.grid(pop1 = factor(pop_levels, pop_levels),
                         pop2 = factor(pop_levels, pop_levels))

lower_triangle_dxy <- full_grid %>%
  mutate(p1_num = as.numeric(pop1), p2_num = as.numeric(pop2)) %>%
  filter(p1_num > p2_num) %>%  # Lower triangle, no diagonal
  left_join(inp.avg, by = c("pop1", "pop2")) %>%
  left_join(inp.avg %>% rename(pop1 = pop2, pop2 = pop1, Mean_rev = Mean), 
            by = c("pop1", "pop2")) %>%
  mutate(Mean = ifelse(is.na(Mean), Mean_rev, Mean)) %>%
  select(-Mean_rev, -p1_num, -p2_num)

# Plot
plot.inp.dxy <- ggplot(lower_triangle_dxy, aes(x = fct_rev(fct_infreq(pop1)), 
                           y = fct_infreq(pop2), fill = Mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Mean)), size = 4) +
  scale_fill_gradient2(low = "grey25", high = "red", mid = "pink", 
                       midpoint = 3.5, limits = c(1.5, 5),
                       name = expression(atop(Avg~dxy, (x10^{-3})))
  ) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_blank()) +
  coord_fixed()

plot.inp.dxy



#####################################################
####################################################
# Read and process FST data
inp.fst <- read.table("AD1380_8groups.fst3.txt", sep = "", header = FALSE) %>% 
  rename(pop1 = V1, pop2 = V2, chromosome = V3, avg_fst = V4, fst_sd = V5) %>%
  mutate(pop1 = recode(pop1, "YUCE" = "YUC-E", "YUCW" = "YUC-W"),
         pop2 = recode(pop2, "YUCE" = "YUC-E", "YUCW" = "YUC-W"))

pop_levels <- c("Cultivar", "LR1", "LR2", "PR", "GD", "FL", "YUC-E", "YUC-W")

# Get averages
inp.fst.avg <- inp.fst %>%
  group_by(pop1, pop2) %>%
  summarise(Mean = mean(avg_fst), .groups = "drop") %>%  # Fst usually not multiplied by 1000
  mutate(pop1 = factor(pop1, pop_levels),
         pop2 = factor(pop2, pop_levels))

# Create full matrix and fill symmetric
full_grid_fst <- expand.grid(pop1 = factor(pop_levels, pop_levels),
                             pop2 = factor(pop_levels, pop_levels))

lower_triangle_fst <- full_grid_fst %>%
  mutate(p1_num = as.numeric(pop1), p2_num = as.numeric(pop2)) %>%
  filter(p1_num > p2_num) %>%  # Lower triangle, no diagonal
  left_join(inp.fst.avg, by = c("pop1", "pop2")) %>%
  left_join(inp.fst.avg %>% rename(pop1 = pop2, pop2 = pop1, Mean_rev = Mean), 
            by = c("pop1", "pop2")) %>%
  mutate(Mean = ifelse(is.na(Mean), Mean_rev, Mean)) %>%
  select(-Mean_rev, -p1_num, -p2_num)

# Plot Fst
plot.inp.fst <- ggplot(lower_triangle_fst, aes(x = fct_rev(fct_infreq(pop1)), 
                                               y = fct_infreq(pop2), fill = Mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Mean, 2)), size = 4) +  # Fst values are smaller, use more decimals
  scale_fill_gradient2(low = "grey25", high = "red", mid = "pink", 
                       midpoint = 0.425,  # Adjust midpoint for Fst range
                       limits = c(0, 0.85),  # Adjust limits for Fst range
                       name = expression(atop(Avg~F[ST], ""))  # Fst doesn't need ×10⁻³
  ) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_blank()) +
  coord_fixed()

plot.inp.fst

############################################################################################
############################################################################################



############################################################################################
############################################################################################

heteperce <- read.csv("AD1_n380.AhDh.combined.bi.het", sep = '\t') %>%
  rename(INDV = INDV, O_HOM = O.HOM., E_HOM = E.HOM., N_SITES = N_SITES, F = F) %>%
  mutate(prefix = str_extract(INDV, "^[^_]+_[^_]+"),
         group = case_when(str_detect(INDV, "^AD1_YUC_RiCh") ~ "AD1_YUC-E",
                           str_detect(INDV, "^AD1_YUC_RiCa") ~ "AD1_YUC-E",
                           str_detect(INDV, "^AD1_YUC_") ~ "AD1_YUC-W",
                           TRUE ~ prefix),
         He = ((N_SITES - O_HOM) / N_SITES)*100)  %>%
  mutate(group = gsub("AD1_","", group)) %>%
  select(-prefix)



heteperce_mean <- heteperce %>%
  group_by(group) %>% 
  dplyr::summarize(percent_mean = mean(He, na.rm=TRUE)) 

group_ordered <- with(heteperce_mean,  reorder(group,  percent_mean,  mean))
heteperce$group <- factor(heteperce$group, levels = levels(group_ordered))

heteperce$group <- factor(heteperce$group,  levels = c("Cultivar", "LR1", "LR2", "GD", "PR","FL", "YUC-E", "YUC-W"))


heteperce_Plot <- ggplot(heteperce, aes(x=group, y=He, fill = group)) + 
  geom_boxplot(show.legend = F, width = 0.7, alpha = 0.5, outlier.shape = NA) + 
#  geom_text(aes(label = formatC(after_stat(y), format = "f", digits  = 1), group = group), size = 3,
#            stat = 'summary', fun = mean,  nudge_y = 3 , color="blue") +
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue", fill="blue") +
  scale_fill_manual(values=safe_colorblind_palette)  +
  geom_jitter(color="black", size=1, alpha=0.9) +
  #xlab("Populations/Groups") +
  ylab("Heterozygous sites (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = -45, v = -0.5, size = 10),
        panel.border = element_rect(colour = "black", fill=NA))
heteperce_Plot


f_plot <- ggplot(heteperce, aes(x = group, y = F, fill = group)) + 
  geom_boxplot(show.legend = F, width = 0.7, alpha = 0.5, outlier.shape = NA) + 
#  geom_text(aes(label = formatC(after_stat(y), format = "f", digits = 3)),  # F values need more decimals
#            stat = 'summary', fun = mean, size = 3, color = "blue",nudge_y = -0.04) +  # Smaller nudge for F
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "blue") +
  scale_fill_manual(values = safe_colorblind_palette) +
  geom_jitter(color = "black", size = 1, alpha = 0.9, width = 0.2) +
  #xlab("Populations/Groups") +
  ylab(expression(Inbreeding~coefficient~(F[IS])))+
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 10),
        panel.border = element_rect(colour = "black", fill = NA))

f_plot
############################################################################################
############################################################################################



############################################################################################
############################################################################################
library(detectRUNS)
load("../Fig3_ROH/AD1_n380.final.RData")

slidingRuns$group <- gsub("AD1_", "", slidingRuns$id)

slidingRuns2 <- slidingRuns %>%
  filter(lengthBps >= 250000) %>%
  mutate(category = case_when(
    lengthBps >= 250000  & lengthBps < 500000  ~ "0.25–0.5 Mb",
    lengthBps >= 500000  & lengthBps < 1000000 ~ "0.5–1 Mb",
    lengthBps >= 1000000 & lengthBps < 2000000 ~ "1–2 Mb",
    lengthBps >= 2000000 ~ "2 Mb +",
    TRUE ~ NA_character_))



slidingRuns4 <- slidingRuns2 %>%
  group_by(group) %>% 
  dplyr::summarize(Froh = sum(lengthBps, na.rm=TRUE)/2296245394) %>%
  mutate(group2 = case_when(
    str_detect(group, "^(YUC_RiCa|YUC_RiCh)") ~ "YUC-E",
    str_detect(group, "^YUC_") ~ "YUC-W",
    TRUE ~ str_extract(group, "^[^_]+"))) %>%  
  as.data.frame()

slidingRuns4$group2 <- factor(slidingRuns4$group2,  levels = c("Cultivar", "LR1", "LR2", "GD", "PR","FL", "YUC-E", "YUC-W"))


mean_Froh <- slidingRuns4 %>%
  group_by(group2) %>%
  summarise(mean_Froh = mean(Froh, na.rm = TRUE)) 

Roh_Froh_Plot <- ggplot(slidingRuns4, aes(x=group2, y=Froh, fill=group2)) + 
  geom_boxplot(show.legend = F, width = 0.7, alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(color="black", size=1,  width = 0.1, alpha=0.7) +
#  geom_text(data = mean_Froh, aes(x = group2, label = formatC(mean_Froh, format = "f", digits = 3)),
#            y = 0.02, size = 4) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue", fill="blue") +
  scale_fill_manual(values=safe_colorblind_palette)+
  xlab("Populations/Groups") +
  ylab(expression("F"[ROH]*" inbreeding")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0, size = 10),
        panel.border = element_rect(colour = "black", fill=NA))

Roh_Froh_Plot


############################################################################################
############################################################################################



############################################################################################
############################################################################################

finalplot <- ggdraw() +
  draw_plot(plot.inp.pi.avg2, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(plot.inp.dxy, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(heteperce_Plot, x = 0, y = 0, width = 0.333333, height = 0.5) +
  draw_plot(f_plot, x = 0.333333, y = 0, width = 0.333333, height = 0.5) +
  draw_plot(Roh_Froh_Plot, x = 0.666666, y = 0, width = 0.333333, height = 0.5) +
  
  draw_plot_label(label = c("a", "b", "c", "d", "e"), size = 17, fontface = "bold",
                  x = c(0, 0.5, 0, 0.333333, 0.666666), y = c(1, 1, 0.5, 0.5, 0.5))

finalplot2 <- ggdraw() +
  draw_plot(plot.inp.pi.avg2, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(plot.inp.dxy, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(heteperce_Plot, x = 0, y = 0, width = 0.333333, height = 0.5) +
  draw_plot(f_plot, x = 0.333333, y = 0, width = 0.333333, height = 0.5) +
  draw_plot(Roh_Froh_Plot, x = 0.666666, y = 0, width = 0.333333, height = 0.5) +
  
  draw_plot_label(label = c("a", "b", "c", "d", "e"), size = 17, fontface = "bold",
                  x = c(0, 0.5, 0, 0.333333, 0.666666), y = c(1, 1, 0.5, 0.5, 0.5))

finalplot2

pdf("../Fig4_PixyAD1380_FIS_HE_ROH.pdf", width = 13, height = 9)
finalplot2
dev.off()

