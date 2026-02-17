setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig4_GeneticLoad/Fig4_sift4g/")


library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

#######################################################################
#######################################################################

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
                        "YUC" = 15,
                        'AD2'= 1,
                        'AD4'= 2)

#######################################################################
#######################################################################

deleterious <- fread("highconf_gt_with_score_outgroup0.tsv")[,-383]

# 查看SIFT分数分布
hist(as.numeric(deleterious$SIFT_SCORE), 
     main = "SIFT Score Distribution",
     xlab = "SIFT Score",
     col = "lightblue")
    
ggplot(deleterious, aes(x = CHROM, y = POS, color = SIFT_SCORE)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_gradient(low = "blue",    # max value 0.04 → red
                       high = "yellow",  # min value 0 → blue    
                       limits = c(0, 0.04),
                       oob = scales::squish,
                       guide = guide_colorbar(reverse = TRUE)  ) +
  theme_classic() +
  labs(  x = "Chromosome",
    y = "Position on Chromosome",
    color = "SIFT score",
    title = "Deleterious Variant Positions Colored by SIFT Score" ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####################################################################

deleterious_long <- deleterious %>%
  pivot_longer(  cols = -c(CHROM, POS, SIFT_SCORE), 
                 names_to = "Sample", 
                 values_to = "Genotype"  ) %>%
  mutate(group = gsub("AD1_", "", Sample)) %>%
  mutate(
    # 危害权重：SIFT越低 → 权重越高
    Harm_Weight = 1 - (as.numeric(SIFT_SCORE) / 0.04),
    
    # Additive model
    Harmful_Alleles_Additive = case_when(
      Genotype == "0/0" ~ 0,
      Genotype == "0/1" | Genotype == "1/0" ~ 0.5,
      Genotype == "1/1" ~ 1,
      TRUE ~ 0),
    # Recessive model
    Harmful_Alleles_Recessive = case_when(
      Genotype == "0/0" ~ 0,
      Genotype == "0/1" | Genotype == "1/0" ~ 0,
      Genotype == "1/1" ~ 1,
      TRUE ~ 0))

# Weighted harmful allele burden
sample_burden_weighted <- deleterious_long %>%
  mutate(weighted_allele_Additive = Harmful_Alleles_Additive * Harm_Weight) %>%  # Multiply by SIFT_SCORE
  mutate(weighted_allele_Recessive = Harmful_Alleles_Recessive * Harm_Weight) %>%  # Multiply by SIFT_SCORE
  #mutate(Subgenome = ifelse(grepl("^A", CHROM), "A", ifelse(grepl("^D", CHROM), "D", NA))) %>% 
  group_by(group) %>%
  summarise( total_weighted_burden_Additive = sum(weighted_allele_Additive, na.rm = TRUE),
             total_weighted_burden_Recessive = sum(weighted_allele_Recessive, na.rm = TRUE),
             .groups = "drop")%>%
  mutate(group2 = if_else(str_detect(group, "^(Cultivar|LR1|LR2|GD)_"),
                          str_extract(group, "^[^_]+"),
                          str_extract(group, "^[^_]+_[^_]+")),
         group3 = case_when(
           str_detect(group, "^(YUC_RiCa|YUC_RiCh)") ~ "YUC-E",
           str_detect(group, "^YUC_") ~ "YUC-W",
           TRUE ~ str_extract(group, "^[^_]+")))  %>%
  pivot_longer( cols = c(total_weighted_burden_Additive, total_weighted_burden_Recessive),
                names_to = "Model",
                values_to = "Weighted_Burden" ) %>%
  mutate(Model = recode(Model,
                        "total_weighted_burden_Additive" = "Additive",
                        "total_weighted_burden_Recessive" = "Recessive"))



group_comparisons <- list(
  c("GD", "YUC-W"),
  c("GD", "YUC-E"),
  c("YUC-W", "YUC-E"))

sample_burden_weighted$group3 <- factor(
  sample_burden_weighted$group3,  
  levels = c("Cultivar", "LR1", "LR2", "GD", "PR","FL", "YUC-E", "YUC-W"))

sift4gplot <- ggplot(sample_burden_weighted,
       aes(x = group3, y = Weighted_Burden, colour = group3, fill = group3)) +
  geom_boxplot(alpha = 0.6, width = 0.9, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8, color = "black") +
  facet_wrap(~ Model, ncol = 1, scales = "free_x", space = "free_y",strip.position = "left") +
  scale_fill_manual(values = safe_colorblind_palette) +
  theme_classic() +
  labs(y = "SIFT4G weighted genetic load (SIFT score < 0.04)") +
  theme( legend.position = "none",
         panel.background = element_blank(),    
         plot.background = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_text(size = 12),
         axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
         axis.text.y = element_text(size = 10),
         strip.placement = "outside",
         strip.text = element_text(size = 10),
         strip.background = element_blank())

##################
##################

# Load necessary library
library(dplyr)

# Make sure group3 is a factor and set YUC-W as the reference
sample_burden_weighted$group3 <- factor(sample_burden_weighted$group3, levels = c("YUC-W",
                                                                "Cultivar", "LR1", "LR2", "GD", "PR", "FL", "YUC-E"))

# Separate the two models
burden_add <- sample_burden_weighted %>% filter(Model == "Additive")
burden_rec <- sample_burden_weighted %>% filter(Model == "Recessive")

# 1. ANOVA for Additive model
anova_add <- aov(Weighted_Burden ~ group3, data = burden_add)
summary(anova_add)

# Tukey post-hoc test for Additive model
tukey_add <- TukeyHSD(anova_add, "group3")
tukey_add

# 2. ANOVA for Recessive model
anova_rec <- aov(Weighted_Burden ~ group3, data = burden_rec)
summary(anova_rec)

# Tukey post-hoc test for Recessive model
tukey_rec <- TukeyHSD(anova_rec, "group3")
tukey_rec
