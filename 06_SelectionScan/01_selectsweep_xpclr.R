setwd("C:/Users/Weixuan/Desktop/YUC_wildAD1/Fig6_WildCultivar_SelectionSweep/Fig6_Part1_selectionsweep/")

library(ggplot2)
library(dplyr)
# Load the library
library(qqman)

# Make the Manhattan plot on the gwasResults dataset
MIN_SNPS <- 10
PI_RATIO_HIGH_Q <- 0.99
PI_RATIO_LOW_Q  <- 0.01
FST_HIGH_Q      <- 0.95

#
#

Cultivar.pi <- read.csv("all_Cultivar_pi.windowed.pi", sep="", header=TRUE)
YUC.pi <- read.csv("all_YUC_pi.windowed.pi", sep="", header=TRUE)

YUC_Cultivar_pi <- Cultivar.pi %>%
  rename(PI_Cultivar=PI, N_VARIANTS_Cultivar=N_VARIANTS) %>%
  inner_join(YUC.pi %>% rename(PI_YUC=PI, N_VARIANTS_YUC=N_VARIANTS), 
             by=c("CHROM","BIN_START","BIN_END")) %>%
  mutate(BIN_MID=(BIN_START+BIN_END)/2, PI_ratio=PI_YUC/PI_Cultivar) %>%
  filter(N_VARIANTS_Cultivar>=10, N_VARIANTS_YUC>=10) %>%
  arrange(CHROM, BIN_START) %>%
  mutate(selection_PI=case_when(
    PI_ratio>=quantile(PI_ratio, 0.95) ~ "Low_diversity_domestication",
    #PI_ratio<=quantile(PI_ratio, 0.05) ~ "High_diversity_introgression",
    TRUE ~ NA_character_))

#
#

YUC_Cultivar_fst <- read.csv("all_fst.windowed.weir.fst", sep="", header=TRUE) %>%
  mutate(BIN_MID=(BIN_START+BIN_END)/2) %>%
  filter(if("N_VARIANTS" %in% colnames(.)) N_VARIANTS>=10 else TRUE) %>%
  mutate(selection_FST=case_when(WEIGHTED_FST>=quantile(WEIGHTED_FST, 0.95) ~ "High_FST", TRUE~NA_character_))

#
#
YUC_Cultivar_xpclr <- read.csv("all_Cultivar_Wild_xpclr.txt", sep="", header=TRUE) %>%
  mutate(BIN_MID=(start+stop)/2) %>%
  filter(if("nSNPs" %in% colnames(.)) nSNPs >=10 else TRUE) %>%
  mutate(selection_xpclr=case_when(xpclr_norm>=quantile(xpclr_norm, 0.95) ~ "High_xpclr", TRUE~NA_character_))

################################################################




ggplot(YUC_Cultivar_pi, aes(x=BIN_MID, y=PI_ratio)) +
  geom_point(size=0.6, alpha=0.7, color="grey50") +
  geom_point(data=YUC_Cultivar_pi %>% filter(!is.na(selection_PI)), aes(color=selection_PI), size=0.8) +
  facet_wrap(~CHROM, ncol=1, scales="fixed", strip.position="right") +
  labs(title="PI Ratio (YUC / Cultivar) with Selection Signals", x="Genomic Position", y="PI Ratio") +
  theme_bw() +
  scale_color_manual(values=c("Low_diversity_domestication"="red","High_diversity_introgression"="blue"))

ggplot(YUC_Cultivar_fst, aes(x=BIN_MID, y=WEIGHTED_FST)) +
  geom_point(size=0.6, alpha=0.7, color="grey50") +
  geom_point(data=YUC_Cultivar_fst %>% filter(!is.na(selection_FST)), aes(color=selection_FST), size=0.8) +
  facet_wrap(~CHROM, ncol=1, scales="fixed", strip.position="right") +
  labs(title="FST (Weighted) with Top 5% Signals", x="Genomic Position", y="Weighted FST") +
  theme_bw() +
  scale_color_manual(values=c("High_FST"="red"))

ggplot(YUC_Cultivar_xpclr, aes(x=BIN_MID, y=xpclr_norm)) +
  geom_point(size=0.6, alpha=0.7, color="grey50") +
  geom_point(data = YUC_Cultivar_xpclr %>% filter(!is.na(selection_xpclr)),
             aes(color=selection_xpclr), size=0.8) +
  facet_wrap(~chrom, ncol=1, scales="fixed", strip.position="right") +
  labs(title="XP-CLR Scores with Selection Signals", x="Genomic Position", y="Normalized XP-CLR") +
  theme_bw() +
  scale_color_manual(values=c("High_xpclr"="red"))
dev.off()

#
#
#
################################################################

combined_select <- YUC_Cultivar_pi %>%
  select(CHROM, BIN_START, BIN_END, BIN_MID, PI_ratio, selection_PI) %>%
  inner_join(
    YUC_Cultivar_fst %>%
      select(CHROM, BIN_START, BIN_END, BIN_MID, WEIGHTED_FST, MEAN_FST, selection_FST),
    by = c("CHROM","BIN_MID"),
    suffix = c("_pi","_fst")) %>%
  left_join(
    YUC_Cultivar_xpclr %>%
      rename(CHROM = chrom) %>%
      select(CHROM, BIN_MID, xpclr_norm, selection_xpclr),
    by = c("CHROM","BIN_MID")  ) %>%
  filter(if_any(c(selection_PI, selection_FST, selection_xpclr), ~ !is.na(.))) %>%
  transmute(
    CHROM, BIN_START = BIN_START_pi, BIN_END = BIN_END_pi, BIN_MID,
    PI_ratio,  WEIGHTED_FST,  MEAN_FST,  xpclr_norm,
    selection_PI,  selection_FST,  selection_xpclr )

combined_select_all_signals <- combined_select %>%
  filter(rowSums(!is.na(select(., selection_PI, selection_FST, selection_xpclr))) >= 2)

nrow(combined_select_all_signals)
colnames(combined_select_all_signals)


combined_select_all_signals <- combined_select_all_signals %>%
  mutate(signal_combination = case_when(
    !is.na(selection_PI) & !is.na(selection_FST) & is.na(selection_xpclr) ~ "PI + FST",
    !is.na(selection_PI) & is.na(selection_FST) & !is.na(selection_xpclr) ~ "PI + XP-CLR",
    is.na(selection_PI) & !is.na(selection_FST) & !is.na(selection_xpclr) ~ "FST + XP-CLR",
    !is.na(selection_PI) & !is.na(selection_FST) & !is.na(selection_xpclr) ~ "PI + FST + XP-CLR")) %>%
  mutate(signal_combination = factor(signal_combination,
                                     levels = c("PI + FST",
                                                "PI + XP-CLR",
                                                "FST + XP-CLR",
                                                "PI + FST + XP-CLR")))
#
#

library(dplyr)
library(ggplot2)
library(tibble)
library(cowplot)
library(stringr)

# 1) Chromosome lengths from faidx
chrom_lengths <- tribble(
  ~CHROM, ~LENGTH,
  "Ah_01",119153486,"Ah_02",108577906,"Ah_03",112634255,"Ah_04",89458773,
  "Ah_05",113043850,"Ah_06",127787289,"Ah_07",100642025,"Ah_08",126877624,
  "Ah_09",85855322,"Ah_10",115166275,"Ah_11",120190126,"Ah_12",107960188,
  "Ah_13",108043650,"Dh_01",66219293,"Dh_02",73967036,"Dh_03",54899334,
  "Dh_04",59297589,"Dh_05",68155215,"Dh_06",67330030,"Dh_07",62911544,
  "Dh_08",72308093,"Dh_09",56184090,"Dh_10",69135794,"Dh_11",74193059,
  "Dh_12",63574961,"Dh_13",65507900)

# 2) Order & map chromosomes to y-axis
chrom_lengths <- chrom_lengths %>%
  mutate(CHROM=str_replace(CHROM,"^Ah_","A"),
         CHROM=str_replace(CHROM,"^Dh_","D"))
chrom_order <- chrom_lengths$CHROM
chrom_lengths <- chrom_lengths %>% mutate(y=match(CHROM,chrom_order))

combined_select_all_signals <- combined_select_all_signals %>%
  mutate(CHROM=str_replace(CHROM,"^Ah_","A"),
         CHROM=str_replace(CHROM,"^Dh_","D"),
         y=match(CHROM,chrom_order))

# 3) Color-blind friendly colors
cb_colors <- c(
  "PI + FST"="#E69F00",
  "PI + XP-CLR"="#56B4E9",
  "FST + XP-CLR"="#009E73",
  "PI + FST + XP-CLR"="#D55E00")

# 4) Plot chromosomes + selection regions
plotdistribution <- ggplot() +
  geom_segment(data = chrom_lengths, aes(x = 0, xend = LENGTH, y = y, yend = y), linewidth = 4, lineend = "butt", color = "grey90", alpha = 0.6) +
  geom_segment(data = combined_select_all_signals, aes(x = BIN_START, xend = BIN_END, y = y, yend = y, color = signal_combination), size = 3, lineend = "butt") +
  scale_color_manual(values = cb_colors) +
  scale_y_continuous(breaks = chrom_lengths$y, labels = chrom_lengths$CHROM, expand = expansion(mult = 0.03)) +
  labs(y = "Chromosome", x = "Genomic position (bp)", color = "Selection signals") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),  # Added angle for y-axis labels
    panel.grid.major.x = element_line(color = "black", linetype = "dashed"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    #axis.title.x = element_text(size = 12, margin = margin(t = 20)),
    #axis.title.y = element_text(size = 12),
    axis.title = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 10),
    legend.position = c(0.77, 0.99),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = NA, color = "black", linewidth = 0.3),
    legend.title = element_text(size = 10, margin = margin(b = 2)),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.margin = margin(3, 3, 3, 3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()  )
  

plotdistribution

######################################################################################



finalplot <- ggdraw() +
  draw_plot(structureplot, x = 0, y = 0, width = 0.4, height = 1) +
  draw_plot(plotdistribution, x = 0.4, y = 0, width = 0.6, height = 1)+ 
  draw_plot_label(label = c("a","b","c"), size = 16, fontface = "bold", x = c(0,0,0.4), y = c(1,0.66,1))
finalplot


pdf("../Fig6_selection_plot.pdf", width = 12, height = 7,    compress = TRUE)
finalplot
dev.off()


png("../Fig6_selection_plot.png", 
    width = 13, 
    height = 7, 
    units = "in",  # 单位设为英寸
    res = 300)     # 分辨率，300dpi适用于大多数出版
finalplot
dev.off()

######################################################################################

#compute length of each window
combined_select_all_signals <- combined_select_all_signals %>%
  mutate(window_length = BIN_END - BIN_START)

# summarize total length per signal_combination
signal_summary <- combined_select_all_signals %>%
  group_by(signal_combination) %>%
  summarise(total_length = sum(window_length)) %>%
  ungroup() %>%
  mutate(fraction = total_length / sum(chrom_lengths$LENGTH))  # compute fraction here

# add total row
signal_summary_total <- bind_rows(
  signal_summary,
  summarise(signal_summary,
            signal_combination = "Total",
            total_length = sum(total_length),
            fraction = sum(fraction))
)

signal_summary_total


# total genome length
genome_length <- sum(chrom_lengths$LENGTH)

# summarize total signal per chromosome and fraction of whole genome
signal_by_chr_any <- combined_select_all_signals %>%
  group_by(CHROM) %>%
  summarise(total_signal = sum(window_length), .groups = "drop") %>%
  mutate(fraction = total_signal / genome_length) %>%
  arrange(desc(total_signal))

signal_by_chr_any
###########################################################################

library(rtracklayer)
genes_gr <- read.delim("AD1.TX2094.v2.gtf", header = FALSE, comment.char = "#", sep = "\t") %>%
  filter(V3 == "gene") %>%
  mutate(V1 = str_replace(V1,"^Ah_","A"),
         V1 = str_replace(V1,"^Dh_","D")) %>%
  { GRanges(seqnames = .$V1,
            ranges  = IRanges(start = .$V4, end = .$V5),
            strand  = .$V7,
            gene_id = .$V9) }

# selection windows as GRanges
selection_gr <- GRanges(
  seqnames = combined_select_all_signals$CHROM,
  ranges   = IRanges(combined_select_all_signals$BIN_START,
                     combined_select_all_signals$BIN_END))

# overlapping genes (unique)
selected_gene_ids <- unique(mcols(genes_gr)$gene_id[subjectHits(findOverlaps(selection_gr, genes_gr))])

# extract proteins
protein_fa <- readAAStringSet("AD1.TX2094.v2.aa.fa")
selected_proteins <- protein_fa[Reduce("|", lapply(selected_gene_ids, \(g) startsWith(names(protein_fa), g)))]

# save
writeXStringSet(selected_proteins, "selected_genes_proteins.fa")


#####################################################################################

selected_windows <- combined_select_all_signals %>% 
  filter(signal_combination == "PI + FST + XP-CLR")

selectionregion <- GRanges(seqnames = selected_windows$CHROM, 
                           ranges = IRanges(start = selected_windows$BIN_START, end = selected_windows$BIN_END))

overlapping_genes_unique <- as.data.frame(genes_gr[subjectHits(findOverlaps(selectionregion, genes_gr))]) %>%
  select(seqnames, start, end, strand, gene_id) %>% 
  distinct(gene_id, .keep_all = TRUE)

# 5️⃣ Extract protein sequences
protein_fa <- readAAStringSet("AD1.TX2094.v2.aa.fa")
idx <- Reduce("|", lapply(overlapping_genes_unique$gene_id, function(g) startsWith(names(protein_fa), g)))
selected_proteins <- protein_fa[idx]

# 6️⃣ Save selected proteins
writeXStringSet(selected_proteins, "selected_genes_PI_FST_XPCLR_proteins.fa")






