setwd("C:/Users/weixuan/Desktop/YUC_wildAD1/Fig4_GeneticLoad/Fig4_SwifTE/")
library(stringr)
library(dplyr)
library(grid)   # for the textGrob() function
library(ggpubr)
library(ggplot2)
library(scales)
library(reshape2)
library(pheatmap)
library("cowplot")
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(readr)

#######################################################################
#######################################################################


######################################################################
### swift ###############################################################
#######################################################################

swift <- read.table("final_n116.bed", sep = "", header = FALSE, comment.char = "") %>% 
  rename(chrom = V1, start = V2, end = V3, 
         TE_name = V4, score = V5, strand = V6, samplename = V7) %>%
  mutate(
    TE_name2 = str_extract(TE_name, "(?<=#)[^,;]+"),
    TE_type = str_split(TE_name2, "/", simplify = TRUE)[, 1],
    TE_type = ifelse(TE_type == "Unknown", "unknown", TE_type),  # Only change Unknown
    TE_family = str_split(TE_name2, "/", simplify = TRUE)[, 2],
    TE_family = ifelse(TE_family == "Unknown", "unknown", TE_family),  # Only change Unknown
    samplename = str_replace(samplename, "^YUC-", "GD_"), 
    samplename = str_remove(samplename, "_merged_R1"),
    Species = case_when(
      str_detect(samplename, "^(YUC_RiCa|YUC_RiCh)") ~ "YUC-E",
      str_detect(samplename, "^YUC_") ~ "YUC-W",
      TRUE ~ str_extract(samplename, "^[^_]+"))) %>%
  arrange(chrom, start)



inds_per_species <- swift %>%
  distinct(Species, samplename) %>%  # or use a column identifying individuals
  count(Species, name = "num_individuals")

# Create summary
swift_summary <- swift %>%
  filter(!TE_type %in% c("unknown", "trna")) %>%  # lowercase now
  
  # Keep positions with ≥2 individuals per species
  add_count(chrom, start, end, TE_type,  name = "n_individuals") %>%
  filter(n_individuals >= 3) %>%
  
  # Count positions per TE type/family
  count(chrom, Species, TE_type, TE_family, name = "Count") %>%
  
  # Add individual counts and calculate percentages
  left_join(inds_per_species, by = "Species") %>%
  mutate( perc = Count / num_individuals,
          Group = str_c(Species, TE_type, sep = "_"),
          chromsub = case_when(
            str_sub(chrom, 1, 1) == "A" ~ "A subgenome",
            str_sub(chrom, 1, 1) == "D" ~ "D subgenome",
            TRUE ~ str_sub(chrom, 1, 1)))  %>%  
  filter(TE_type %in% c("LTR", "DNA")) 


swift_summary$TE_type <- factor(swift_summary$TE_type, levels =  c("LTR","LINE","SINE", "DNA","MITE"))
swift_summary$TE_family <- factor(swift_summary$TE_family, levels=c("Copia","Gypsy","unknown","DTA","DTC","DTH","DTM","DTT","MULE-MuDR","CMC-EnSpm","hAT-Tip100","PIF-Harbinger","hAT-Tag1","Helitron"))

colors18 <- c(
  brewer.pal(4, "Set1"),
  brewer.pal(4, "Set2"),
  brewer.pal(6, "Set3") )

swift_plot_percent <- ggplot(swift_summary, aes(x = Species, y = perc, fill = TE_family)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = colors18) +
  labs(fill = "TE Superfamily", y = "TIP count per individual") +
  facet_grid(TE_type ~ chromsub, scales = "free_y", switch = "both") +
  guides(fill = guide_legend(ncol = 2, byrow = FALSE)) +  # 2-column legend
  theme(
    panel.background = element_blank(),    
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.placement = "outside",
    strip.text = element_text(size = 10),
    # Position legend inside plot
    legend.position = c(0.47, 0.99),  # Top-left corner
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA, color = "black", linewidth = 0.3),
    legend.title = element_text(size = 8, margin = margin(b = 2)),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.05, "cm"),  # Less vertical space
    legend.spacing.x = unit(0.1, "cm"),   # Space between columns
    legend.margin = margin(3, 3, 3, 3))


print(swift_plot_percent)


pdf("Figxx_swift_tips2_cotton.pdf", width = 5, height = 12)
swift_plot_percent
dev.off()


######################################################################################################
######################################################################################################
######################################################################################################

ltr <- swift_summary %>%
  group_by(chromsub, TE_type, Group) %>%
  summarise(Total_Count = sum(perc), .groups = "drop") %>%
  filter(TE_type == "LTR")

# pivot wider to compare A vs D
A_vs_D <- ltr %>%
  select(chromsub, Group, Total_Count) %>%
  pivot_wider(
    names_from = chromsub,
    values_from = Total_Count) %>%
  mutate(
    diff = `A-genome` - `D-genome`,
    ratio = `A-genome` / `D-genome`)

A_vs_D

# A-genome vs YUC-W
A_ref <- ltr %>% filter(chromsub == "A-genome", Group == "YUC-W_LTR") %>% pull(Total_Count)

A_compare <- ltr %>%
  filter(chromsub == "A-genome") %>%
  mutate(
    fold_change_vs_YUCW = Total_Count / A_ref,
    diff_vs_YUCW = Total_Count - A_ref
  )

# D-genome vs YUC-W
D_ref <- ltr %>% filter(chromsub == "D-genome", Group == "YUC-W_LTR") %>% pull(Total_Count)

D_compare <- ltr %>%
  filter(chromsub == "D-genome") %>%
  mutate(
    fold_change_vs_YUCW = Total_Count / D_ref,
    diff_vs_YUCW = Total_Count - D_ref
  )

A_compare
D_compare


#



















library(UpSetR)
library(ggplot2)

library(UpSetR)
library(ggplot2)

# Prepare data for upset plot
upset_data <- swift %>%
  filter(!TE_type %in% c("unknown", "trna")) %>%
  add_count(chrom, start, end, TE_type, Species, name = "n_individuals") %>%
  filter(n_individuals >= 2) %>%
  
  # Create binary matrix: 1 if species has TE at position
  distinct(chrom, start, end, TE_type, Species) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols = c(chrom, start, end, TE_type),
    names_from = Species,
    values_from = value,
    values_fill = 0
  )

# Convert to matrix format for UpSetR
upset_matrix <- as.data.frame(upset_data[, -c(1:3)])

# Create upset plot
upset(upset_matrix, 
      nsets = ncol(upset_matrix),
      nintersects = 20,
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "darkred")


#






##### merge overlapped genomic regions in TEs
library(GenomicRanges)

count_TE_gr <- GRanges(seqnames = count_TE$chrom,
              ranges = IRanges(start = count_TE$start, end = count_TE$end),
              yuc = count_TE$yuc_non_na,
              fl = count_TE$fl_non_na,
              pr = count_TE$pr_non_na,
              gd = count_TE$gd_non_na)

count_TE_reduced_gr <- reduce(count_TE_gr, with.revmap = TRUE)

sums <- lapply(mcols(count_TE_reduced_gr)$revmap, function(idx) {
  colSums(as.data.frame(mcols(count_TE_gr)[idx, c("yuc", "fl", "pr", "gd")]))
})

sum_df <- do.call(rbind, sums)

mcols(count_TE_reduced_gr) <- cbind(mcols(count_TE_reduced_gr), sum_df)



















###################### finding relationships to nearest gene
#BiocManager::install("rtracklayer")

library(rtracklayer)

gtf_raw <- read.delim("AD1.TX2094.v2.gtf", header = FALSE, comment.char = "#", sep = "\t")
genes_df <- gtf_raw[gtf_raw$V3 == "gene", ]

genes_gr <- GRanges(
  seqnames = genes_df$V1,
  ranges = IRanges(start = genes_df$V4, end = genes_df$V5),
  strand = genes_df$V7,
  gene_id = genes_df$V9  # raw gene name from column 9
)

hits <- distanceToNearest(count_TE_reduced_gr, genes_gr)

count_TE_reduced_gr$nearest_gene <- mcols(genes_gr)[subjectHits(hits), "gene_id"]
count_TE_reduced_gr$distance_to_gene <- mcols(hits)$distance

count_TE_reduced_gr$distance_category <- cut(
  count_TE_reduced_gr$distance_to_gene,
  breaks = c(-1, 0, 1000, 2000, 5000, Inf),
  labels = c("within genic", "<1 kbp", "1–2 kbp", "2-5 kbp", ">5 kbp"),
  right = TRUE)

# Compute the sum of the four columns for each row
count_TE_reduced_gr$total_counts <- rowSums(as.data.frame(mcols(count_TE_reduced_gr)[, c("yuc", "fl", "pr", "gd")]))

# Filter GRanges object to keep only rows with total sum >= 100
count_TE_reduced_gr2 <- count_TE_reduced_gr[count_TE_reduced_gr$total_counts >= 20]

count_TE_reduced_gr2[count_TE_reduced_gr2$distance_category == "within genic"]

table(count_TE_reduced_gr2$distance_category)

mcols(count_TE_reduced_gr2)












