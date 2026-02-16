setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig3_Pixy_He_FIS_ROH/Fig3_ROH/")

library(detectRUNS)
library(ggplot2)
library(tidyverse)
library(reshape2)

safe_colorblind_palette <- c("LR1" = "#F4A582",
                             "LR2" = "#D6604D",
                             "Cultivar" = "#B2182B",
                             "PR"= "black",
                             "GD" = "#4D9221",
                             "FL" = "#3B8ABE",
                             "YUC-E" = "#8C66AF",
                             "YUC-W" = "#8C66AF")
safe_shape_palette <- c("FL"= 19,
                        "LR1" = 17,
                        "Cultivar" = 17,
                        "LR2" = 17,
                        "GD" = 3,
                        "PR" = 13,
                        "YUC-W" = 15,
                        "YUC-E" = 1)

#################################################################################
#################################################################################

load("AD1_n380.final.RData")

slidingRuns$group <- gsub("AD1_", "", slidingRuns$id)

slidingRuns2 <- slidingRuns %>%
  filter(lengthBps >= 250000) %>%
  mutate(category = case_when(
    lengthBps >= 250000  & lengthBps < 500000  ~ "0.25–0.5 Mb",
    lengthBps >= 500000  & lengthBps < 1000000 ~ "0.5–1 Mb",
    lengthBps >= 1000000 & lengthBps < 2000000 ~ "1–2 Mb",
    lengthBps >= 2000000 ~ "2 Mb +",
    TRUE ~ NA_character_))


slidingRuns3 <- slidingRuns2 %>%
  group_by(category, group) %>%
  summarise(lengthBps_sum = sum(lengthBps, na.rm = TRUE) / 2296245394, .groups = "drop") %>%
  mutate(group2 = if_else(str_detect(group, "^(Cultivar|LR1|LR2|GD)_"),
                     str_extract(group, "^[^_]+"),
                     str_extract(group, "^[^_]+_[^_]+")),
         group3 = case_when(
           str_detect(group, "^(YUC_RiCa|YUC_RiCh)") ~ "YUC-E",
           str_detect(group, "^YUC_") ~ "YUC-W",
           TRUE ~ str_extract(group, "^[^_]+"))) %>%
  as.data.frame()


category_fac <- as.factor(slidingRuns3$category)
slidingRuns3$category <- factor(slidingRuns3$category, levels = c("2 Mb +","1–2 Mb","0.5–1 Mb","0.25–0.5 Mb"))
slidingRuns3$group2 <- factor(slidingRuns3$group2,  
                              levels = c("Cultivar","LR2","LR1","GD",
                                         "PR_Ph",
                                         "FL_CPH","FL_CPT","FL_FMY","FL_MK","FL_NP","FL_OPB2","FL_OPB4","FL_OPB5","FL_PK","FL_RBD","FL_RBT","FL_RNRB","FL_SR","FL_TC","FL_VKPA",
                                         "YUC_CeCo","YUC_CeDo","YUC_CeDr","YUC_CePr","YUC_RiCa","YUC_RiCh","YUC_SiPa","YUC_SiPr"))

ROH_plot <- ggplot(slidingRuns3,  aes(fill = category, x = group, y = lengthBps_sum)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = c( "0.25–0.5 Mb" = "steelblue",
                "0.5–1 Mb"  = "#56B4E9",
                "1–2 Mb"  = "#E69F00",
                "2 Mb +"  = "#009E73"), drop = FALSE) +
  ylab("Proportion of total genome") +
  facet_wrap(~group2, scales = "free_x", space = "free_x", nrow = 1, strip.position = "bottom") +
  theme_classic() +
  theme(legend.position = c(0.95, 0.85),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),   # hide x-axis labels
    axis.ticks.x = element_blank(),  # hide x-axis ticks
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", colour = "grey30"),
    strip.background = element_blank(),  # remove facet box
    strip.text = element_text(size = 8,  angle = -90)) +
  theme(strip.placement = "outside")


ROH_plot$labels$fill <- "Length of ROH"

ROH_plot

#####################################################################
#####################################################################


slidingRuns3 <- slidingRuns2 %>%
  group_by(category, group) %>%
  summarise(lengthBps_sum = sum(lengthBps, na.rm = TRUE) / 2296245394, .groups = "drop") %>%
  mutate(group2 = if_else(str_detect(group, "^(Cultivar|LR1|LR2|GD)_"),
                          str_extract(group, "^[^_]+"),
                          str_extract(group, "^[^_]+_[^_]+")),
         group3 = case_when(
           str_detect(group, "^(YUC_RiCa|YUC_RiCh)") ~ "YUC-E",
           str_detect(group, "^YUC_") ~ "YUC-W",
           TRUE ~ str_extract(group, "^[^_]+"))) %>%
  as.data.frame()

ROH_total <- slidingRuns3 %>%
  group_by(group3, category) %>% 
  dplyr::summarize(ROH_sum = mean(lengthBps_sum, na.rm=TRUE)) %>%
  as.data.frame()

group_ordered <- with(ROH_total,  reorder(group3,  ROH_sum,  sum, decreasing = TRUE))
ROH_total$group3 <- factor(ROH_total$group3, levels = levels(group_ordered))
ROH_total$group3 <- factor(ROH_total$group3,  levels = c("Cultivar", "LR1", "LR2", "GD", "PR","FL", "YUC-E", "YUC-W"))

ROH_totals <- ROH_total %>%
  group_by(group3) %>%
  summarise(total_ROH = sum(ROH_sum, na.rm = TRUE)) %>%
  ungroup()

ROH_total$category <- factor(ROH_total$category, levels = c("2 Mb +","1–2 Mb","0.5–1 Mb","0.25–0.5 Mb"))

ROH_plot2 <- ggplot(ROH_total, aes(fill = category, y = group3, x = ROH_sum)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = c( "0.25–0.5 Mb" = "steelblue",
                "0.5–1 Mb"  = "#56B4E9",
                "1–2 Mb"  = "#E69F00",
                "2 Mb +"  = "#009E73"), drop = FALSE) +
  geom_text(data = ROH_totals, aes(x = total_ROH, y = group3, label = round(total_ROH, 3)), hjust = -0.1, size = 4, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0, 0.18)) +
  ylab("Populations/Groups") + xlab("Proportion of total genome") +
  theme_classic() +
  theme(legend.position = c(0.86, 0.85), 
        axis.text = element_text(size = 10), 
        legend.background = element_blank(),   # removes legend background
        panel.background = element_blank(),    # removes panel background
        plot.background = element_blank(),     # removes overall plot background        axis.title.y = element_blank(), 
        legend.key = element_rect(fill = "white", colour = "grey30"))

ROH_plot2$labels$fill <- "Length of ROH"
ROH_plot2



#####################################################################
#####################################################################


slidingRuns5 <- slidingRuns2 %>%
  group_by(group) %>% 
  add_count(group) %>%
  group_by(group,n) %>% 
  dplyr::summarize(ROH_sum = sum(lengthBps, na.rm=TRUE)/2296245394) %>%
  mutate(group2 = case_when(
    str_detect(group, "^(YUC_RiCa|YUC_RiCh)") ~ "YUC-E",
    str_detect(group, "^YUC_") ~ "YUC-W",
    TRUE ~ str_extract(group, "^[^_]+"))) %>%  
  as.data.frame()

slidingRuns5 %>%
  group_by(group2) %>% 
  dplyr::summarize(ROH_sum_mean = mean(ROH_sum, na.rm=TRUE), 
                   n_mean = mean(n, na.rm=TRUE))

cor(slidingRuns5$n, slidingRuns5$ROH_sum)


NR_ROH <- ggplot(slidingRuns5, aes(n, ROH_sum, col = group2, shape = group2)) + 
  geom_point(size=2) +
  xlab("Number of ROH") +
  ylab("Proportion of total genome") + 
  scale_color_manual(values= safe_colorblind_palette) +
  scale_shape_manual(values=safe_shape_palette) +
  theme_bw() +  
  theme(legend.position =  c(0.85, 0.35),
        legend.background = element_blank(),
        legend.key = element_rect(fill = "white", colour = "grey30", linetype="dotted")  ) 

NR_ROH$labels$colour <- "Population/Group"
NR_ROH$labels$shape <- "Population/Group"


NR_ROH


#################################################################################
#################################################################################

library(cowplot)

finalplot <- ggdraw() +
  draw_plot(ROH_plot, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(ROH_plot2, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(NR_ROH, x = 0.5, y = 0, width = 0.5, height = 0.5) +

  
  draw_plot_label(label = c("a", "b", "c"), size = 17, fontface = "bold",
                  x = c(0, 0, 0.5), y = c(1, 0.5, 0.5))

finalplot

pdf("../FigS6_ROH_AD1_n380output.pdf", width = 14, height = 10)
finalplot
dev.off()

#####################################################################
#####################################################################






