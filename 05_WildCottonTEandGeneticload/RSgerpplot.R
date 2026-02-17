setwd("C:/Users/Weixuan/Desktop/YUC_wildAD1/Fig4_GeneticLoad/Fig4_gerp/")
#install.packages(c("data.table", "ggplot2", "stringr"))
#install.packages("RColorBrewer")  # if not installed

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

######################################################################
######################################################################

library(data.table)
library(ggsignif)

# 1. Read genotype matrix
geno_file <- "08_n381_merged_GT_gerp_outgroupfilter.tsv"
geno <- fread(geno_file, header = TRUE)[,-c(1:4,388,389)] %>%
  pivot_longer(  cols = -c(CHROM, POS, Gerp), 
                 names_to = "Sample", 
                 values_to = "Genotype"  ) %>%
  mutate(group = gsub("AD1_", "", Sample)) %>%
  mutate(
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


sample_burden <- geno %>%
  # add Subgenome here (or do this earlier in your pipeline)
  mutate(
    Subgenome = case_when(
      grepl("^Ah_", CHROM) ~ "A-genome",
      grepl("^Dh_", CHROM) ~ "D-genome",
      TRUE ~ "Unknown"),
    allele_Additive  = Harmful_Alleles_Additive  * Gerp,
    allele_Recessive = Harmful_Alleles_Recessive * Gerp) %>%
  
  # <-- include Subgenome in the grouping
  group_by(group) %>%
  summarise(
    total_burden_Additive  = sum(allele_Additive,  na.rm = TRUE),
    total_burden_Recessive = sum(allele_Recessive, na.rm = TRUE),
    .groups = "drop") %>%
  
  mutate(
    group2 = if_else(
      str_detect(group, "^(Cultivar|LR1|LR2|GD)_"),
      str_extract(group, "^[^_]+"),
      str_extract(group, "^[^_]+_[^_]+")),
    group3 = case_when(
      str_detect(group, "^(YUC_RiCa|YUC_RiCh)") ~ "YUC-E",
      str_detect(group, "^YUC_") ~ "YUC-W",
      TRUE ~ str_extract(group, "^[^_]+"))) %>%
  
  pivot_longer( cols = c(total_burden_Additive, total_burden_Recessive),
                names_to = "Model", values_to = "Burden") %>%
  mutate( Model = recode(Model, "total_burden_Additive"  = "Additive",
                         "total_burden_Recessive" = "Recessive"))

rm(geno)

unique(sample_burden$group3)

sample_burden$group3 <- factor(
  sample_burden$group3,  
  levels = c("Cultivar", "LR1", "LR2", "GD", "PR","FL", "YUC-E", "YUC-W"))

Geprplot <- ggplot(sample_burden,
                     aes(x = group3, y = Burden, colour = group3, fill = group3)) +
  geom_boxplot(alpha = 0.6, width = 0.9, show.legend = FALSE, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.8, color = "black") +
  #facet_wrap(Subgenome ~ Model, strip.position = "left") +
  facet_wrap(~ Model, ncol = 1, scales = "free_x", space = "free_y",strip.position = "left") +
  
  scale_fill_manual(values = safe_colorblind_palette) +
  theme_classic() +
  labs(y = "GERP++ genetic load (GERP score > 4)") +
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


Geprplot

###################
##################


# Load necessary library
library(dplyr)

# Make sure group3 is a factor and set YUC-W as the reference
sample_burden$group3 <- factor(sample_burden$group3, levels = c("YUC-W",
                                                                "Cultivar", "LR1", "LR2", "GD", "PR", "FL", "YUC-E"))

# Separate the two models
burden_add <- sample_burden %>% filter(Model == "Additive")
burden_rec <- sample_burden %>% filter(Model == "Recessive")

# 1. ANOVA for Additive model
anova_add <- aov(Burden ~ group3, data = burden_add)
summary(anova_add)

# Tukey post-hoc test for Additive model
tukey_add <- TukeyHSD(anova_add, "group3")
tukey_add

# 2. ANOVA for Recessive model
anova_rec <- aov(Burden ~ group3, data = burden_rec)
summary(anova_rec)

# Tukey post-hoc test for Recessive model
tukey_rec <- TukeyHSD(anova_rec, "group3")
tukey_rec
