setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig4_GeneticLoad/Fig4_00_positions")

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggVennDiagram)
library(eulerr)
library(gridExtra)

################################################################################

#cut -f1,2 highconf_gt_with_score.tsv  > highconf_gt_with_score_positionONLY.tsv
#cut -f1,2,3 04_merged_GT_gerp.tsv > 05_GERP_VCFpos_4plus_withinVCF.txt
#grep -v "#" high_impact.vcf | cut -f1,2 > high_impact_positionsONLY.txt

gerppos <- read.csv("09_GERP_VCFpos_4plus_withinVCF_outgroup0.txt", sep = "") %>%
  setNames(c("Chrom","Pos")) %>%
  transmute(Chrom, Pos, Type1 = "Gerp")

snpeffpos <- read.csv("snpeff_high_impact_n381_outgroup0_positiononly.txt", sep = "", header = F) %>%
  setNames(c("Chrom", "Pos"))%>%
  transmute(Chrom, Pos, Type2 = "SnpEFF")

sift4gpos <- read.csv("sift4g_highconf_gt_with_score_outgroup0_positionsONLY.txt", sep="")%>%
  setNames(c("Chrom", "Pos"))%>%
  transmute(Chrom, Pos, Type3 = "SIFT4G")

merged_positions <- gerppos %>%
  full_join(snpeffpos, by = c("Chrom", "Pos")) %>%
  full_join(sift4gpos, by = c("Chrom", "Pos"))%>%
  mutate(
    Subgenome = case_when(
      grepl("^Ah_", Chrom) ~ "A-genome",
      grepl("^Dh_", Chrom) ~ "D-genome",
      TRUE ~ "Unknown"),
    # Create binary columns
    "GERP++" = ifelse(Type1 == "GERP++", 1, 0),
    SnpEFF = ifelse(Type2 == "SnpEFF", 1, 0),
    SIFT4G = ifelse(Type3 == "SIFT4G", 1, 0)  ) %>%
  mutate(Site = paste(Chrom, Pos, sep = "_"))

# Split sets by subgenome
A_sets <- list(
  "GERP++"   = merged_positions$Site[merged_positions$Subgenome=="A-genome" & !is.na(merged_positions$Type1)],
  SnpEFF = merged_positions$Site[merged_positions$Subgenome=="A-genome" & !is.na(merged_positions$Type2)],
  SIFT4G = merged_positions$Site[merged_positions$Subgenome=="A-genome" & !is.na(merged_positions$Type3)]
)

D_sets <- list(
  "GERP++"   = merged_positions$Site[merged_positions$Subgenome=="D-genome" & !is.na(merged_positions$Type1)],
  SnpEFF = merged_positions$Site[merged_positions$Subgenome=="D-genome" & !is.na(merged_positions$Type2)],
  SIFT4G = merged_positions$Site[merged_positions$Subgenome=="D-genome" & !is.na(merged_positions$Type3)]
)

# Prepare plots
pA <- euler(A_sets)
pD <- euler(D_sets)

# Function to generate labels with names + counts
make_labels <- function(euler_obj) {
  # original.values is a named vector with counts for all regions (single and overlaps)
  counts <- round(euler_obj$original.values)
  # names(euler_obj$original.values) gives the combination names, e.g., "GERP", "GERP&SIFT4G", etc.
  labels <- paste0(names(counts), "\n", counts)
  return(labels)
}

# Plot stacked, showing all counts
final_ven <- grid.arrange(
  plot(pA, fills = c("#3182bd80","#31a35480","#756bb180"), labels = make_labels(pA),main = "A subgenome", main.cex = 1),
  plot(pD, fills = c("#3182bd80","#31a35480","#756bb180"), labels = make_labels(pD), main = "D subgenome", main.cex = 1),
  ncol = 2)

################################################################################


finalplot <- ggdraw() +
  draw_plot(swift_plot_percent, x = 0, y = 0, width = 1/3, height = 1) +
  draw_plot(final_ven, x = 1/3, y = 2/3, width = 2/3, height = 1/3) + 
  draw_plot(sift4gplot, x = 1/3, y = 0, width = 1/3, height = 2/3) +
  draw_plot(Geprplot, x = 2/3, y = 0, width = 1/3, height = 2/3) +
  #draw_plot(ROH_plot2, x = 0.666666, y = 0, width = 0.333333, height = 0.5) +
  
  draw_plot_label(label = c("a", "b", "c", "d"), size = 16, fontface = "bold", 
                  x= c(0, 0.333333333, 1/3, 2/3), y = c(1, 1,2/3, 2/3))

finalplot



pdf("../Fig5_TE_Gerp_SnpEFF_SIFT4G.pdf", width = 13, height = 9)
finalplot
dev.off()

save.image()
