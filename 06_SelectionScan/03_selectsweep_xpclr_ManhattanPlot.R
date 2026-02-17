library(dplyr); library(ggplot2); library(patchwork)

# -----------------------
# 1️⃣ Prepare PI-ratio data
df_pi <- YUC_Cultivar_pi %>%
  mutate(CHROM = gsub("h_", "", CHROM)) %>%   # <-- remove prefix here
  arrange(CHROM, BIN_MID) %>%
  group_by(CHROM) %>%
  mutate(chr_len = max(BIN_MID)) %>%
  ungroup() %>%
  mutate(tot = cumsum(lag(chr_len, default=0)), bp_cum = BIN_MID + tot)

axis_set <- df_pi %>% group_by(CHROM) %>% summarise(center = mean(bp_cum), end = max(bp_cum)) %>% mutate(chr = CHROM)
threshold_pi <- quantile(df_pi$PI_ratio, 0.95)
chr_colors <- rep(c("black","grey30"), length.out = length(unique(df_pi$CHROM)))
names(chr_colors) <- unique(df_pi$CHROM)

manh_pi <- ggplot(df_pi, aes(x = bp_cum, y = PI_ratio, color = CHROM)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_hline(yintercept = threshold_pi, linetype = "dashed", color = "red") +
  geom_vline(xintercept = axis_set$end[-nrow(axis_set)], linetype="dashed", color="grey70") +
  scale_color_manual(values = chr_colors) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  labs(y = "PI ratio") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8))

# -----------------------
# 2️⃣ Prepare FST data
df_fst <- YUC_Cultivar_fst %>%
  mutate(CHROM = gsub("h_", "", CHROM)) %>%   # <-- remove prefix here
  arrange(CHROM, BIN_MID) %>%
  group_by(CHROM) %>%
  mutate(chr_len = max(BIN_MID)) %>%
  ungroup() %>%
  mutate(tot = cumsum(lag(chr_len, default=0)), bp_cum = BIN_MID + tot)

threshold_fst <- quantile(df_fst$WEIGHTED_FST, 0.95)
manh_fst <- ggplot(df_fst, aes(x = bp_cum, y = WEIGHTED_FST, color = CHROM)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_hline(yintercept = threshold_fst, linetype = "dashed", color = "red") +
  geom_vline(xintercept = axis_set$end[-nrow(axis_set)], linetype="dashed", color="grey70") +
  scale_color_manual(values = chr_colors) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  labs(y = "FST") +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8))

# -----------------------
# 3️⃣ Prepare XP-CLR data
df_xpclr <- YUC_Cultivar_xpclr %>%
  mutate(chrom = gsub("h_", "", chrom)) %>%   # <-- remove prefix here
  arrange(chrom, BIN_MID) %>%
  group_by(chrom) %>%
  mutate(chr_len = max(BIN_MID)) %>%
  ungroup() %>%
  mutate(tot = cumsum(lag(chr_len, default=0)), bp_cum = BIN_MID + tot)

threshold_xpclr <- quantile(df_xpclr$xpclr_norm, 0.95)
manh_xpclr <- ggplot(df_xpclr, aes(x = bp_cum, y = xpclr_norm, color = chrom)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_hline(yintercept = threshold_xpclr, linetype = "dashed", color = "red") +
  geom_vline(xintercept = axis_set$end[-nrow(axis_set)], linetype="dashed", color="grey70") +
  scale_color_manual(values = chr_colors) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  labs(x = "Chromosome", y = "XP-CLR",) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle=60, hjust=1, size=8),
        axis.text.y = element_text(size=8))

# -----------------------
# Combine plots in one column
manh_all <- manh_pi / manh_fst / manh_xpclr
manh_all


pdf("../FigS10_Manh_selection_plot.pdf", width = 13, height = 7,    compress = TRUE)
manh_all
dev.off()
