library(dbplyr)
library(stringr)

n158kmer <- read.csv("n158.intersections.rename.txt", sep = "", header =T)

sample_cols <- colnames(n158kmer)[-(1:2)]
group_names <- str_extract(sample_cols, "^([^_]+_[^_]+)")
group_map <- split(sample_cols, group_names)

group_means <- sapply(group_map, function(cols) {
  rowMeans(n158kmer[, cols], na.rm = TRUE)})

# Define groups in the desired order
groups_order <- c("AD1_Cultivar", "AD1_LR1", "AD1_LR2", "AD1_YUC", "AD1_FL",
                  "AD1_GD", "AD1_PR", "AD2_Wild", "AD4_mus")

result <- cbind(n158kmer[,1:2], as.data.frame(group_means))%>% 
  filter(!grepl("AD", .[[2]])) %>%
  mutate(across(all_of(groups_order), ~ ifelse(. < 1, 0, .)))


#################################################

# Define AD1, AD2, AD4 groups
ad1_groups <- c("AD1_Cultivar", "AD1_LR1", "AD1_LR2", "AD1_YUC", "AD1_GD", "AD1_PR", "AD1_FL")
ad1_focus <- c("AD1_YUC", "AD1_GD", "AD1_PR", "AD1_FL")
ad2_group <- "AD2_Wild"
ad4_group <- "AD4_mus"
other_groups <- setdiff(groups_order, ad1_focus)


result_long <- result %>%
  # Remove rows where all groups are > 0
  filter(rowSums(across(all_of(groups_order), ~ . > 0)) < length(groups_order)) %>%
  
  # Remove rows where (any AD1>0 AND AD4>0 AND AD2==0)
  filter(!((rowSums(across(all_of(ad1_groups), ~ . > 0)) > 0) & 
             .data[[ad4_group]] > 0 & 
             .data[[ad2_group]] == 0)) %>%
  
  # Remove rows where (all AD1>0 AND AD2==0 AND AD4==0)
  filter(!((rowSums(across(all_of(ad1_groups), ~ . > 0)) == length(ad1_groups)) &
             .data[[ad2_group]] == 0 & 
             .data[[ad4_group]] == 0)) %>%
  
  # Remove rows where (AD2>0 AND all AD1==0 AND AD4==0)
  filter(!( .data[[ad2_group]] > 0 &
              rowSums(across(all_of(ad1_groups), ~ . > 0)) == 0 &
              .data[[ad4_group]] == 0)) %>%
  
  # Remove rows where (AD4>0 AND AD2==0 AND all AD1==0)
  filter(!( .data[[ad4_group]] > 0 &
              .data[[ad2_group]] == 0 &
              rowSums(across(all_of(ad1_groups), ~ . > 0)) == 0)) %>%
  
  # New filter: Remove rows where (AD4==0 AND all AD1>0 AND AD2>0)
  filter(!( .data[[ad4_group]] == 0 &
              rowSums(across(all_of(ad1_groups), ~ . > 0)) == length(ad1_groups) &
              .data[[ad2_group]] > 0)) %>%
  
  # Remove rows where AD2>0 AND AD4>0 AND all AD1==0
  filter(!( .data[[ad2_group]] > 0 &
              .data[[ad4_group]] > 0 &
              rowSums(across(all_of(ad1_groups), ~ . > 0)) == 0)) %>%
  
  # Remove rows where all groups are 0
  filter(rowSums(across(all_of(groups_order), ~ . != 0)) > 0) %>%
  
  # Remove rows where only one group is > 0 (i.e., all other groups are 0)
  # Keep rows where either:
  # 1) more than one group is > 0
  # 2) exactly one AD1-focus group > 10
  filter(
    rowSums(across(all_of(groups_order), ~ . > 0)) > 1 |
      (rowSums(across(all_of(ad1_focus), ~ . > 10)) == 1 &
         rowSums(across(all_of(other_groups), ~ . != 0)) == 0)) %>%
  
  
  # New filter: Remove rows where mean across AD1+AD2+AD4 < 10
  # Calculate mean of non-zero values only
  filter({
    vals <- as.matrix(across(all_of(groups_order)))      # convert to matrix
    rowMeans(ifelse(vals == 0, NA, vals), na.rm = TRUE) >= 10  }) %>%
  
  # Continue with pivot and category assignment
  mutate(Sample = paste0(Number_of_kmers, "_", fasta_File)) %>%
  select(Sample, all_of(groups_order)) %>%
  pivot_longer(cols = all_of(groups_order), names_to = "Group", values_to = "Value") %>%
  mutate(Category = case_when(
    Value == 0 ~ "0",            
    Value < 10 ~ "<10",
    Value >= 10 & Value < 20 ~ "10-20",
    Value >= 20 & Value < 30 ~ "20-30",
    Value >= 30 ~ "30+"
  ),
  Category = factor(Category, levels = c("0", "<10", "10-20", "20-30", "30+")),
  Group = factor(Group, levels = rev(groups_order)))

#########################################################
cb_palette <- c("<10" = "#E69F00",   # orange
                "10-20" = "#56B4E9", # blue
                "20-30" = "#009E73", # green
                "30+" = "#D55E00")   # red

#########################################################
# Compute overlap pattern per sample
sample_order <- result_long %>%
  group_by(Sample) %>%
  summarise(
    n_groups = sum(Value > 0),                             # number of groups >0
    pattern = paste(Group[Value > 0], collapse = "_")      # concatenated non-zero groups
  ) %>%
  arrange(n_groups, pattern) %>%                            # sort by # of groups, then pattern
  pull(Sample)

#######################################

desired_order <- c(
  "87_i366.fasta",   "101_i296.fasta", "459_i56.fasta", "97_i306.fasta", "167_i178.fasta", "82_i396.fasta",
  "484_i55.fasta","193_i153.fasta", "319_i81.fasta", "107_i277.fasta", "116_i255.fasta", "131_i229.fasta",
  "89_i354.fasta", "101_i297.fasta", 
  "417_i64.fasta", "151_i203.fasta", "193_i152.fasta", "149_i204.fasta", "4643_i13.fasta",
  "211_i137.fasta","358_i70.fasta","104_i283.fasta" )

# Apply this order to the Sample factor
result_long$Sample <- factor(result_long$Sample, levels = desired_order)


########################################

sample_groups <- c(
  "87_i366.fasta" = 1, "101_i296.fasta" = 1, "459_i56.fasta" = 1, "97_i306.fasta" = 1, "167_i178.fasta" = 1, "82_i396.fasta" = 1,
  "484_i55.fasta" = 2, "193_i153.fasta" = 2, "319_i81.fasta" = 2, "107_i277.fasta" = 2, "116_i255.fasta" = 2, "131_i229.fasta" = 2,
  "89_i354.fasta" = 3, "101_i297.fasta" = 3,
  "417_i64.fasta" = 4, "151_i203.fasta" = 4, "193_i152.fasta" = 4, "149_i204.fasta" = 4, "4643_i13.fasta" = 4,
  "211_i137.fasta" = 4, "358_i70.fasta" = 4, "104_i283.fasta" = 4)

sample_groups <- sapply(sample_groups, function(x) paste0("Set", x))

# Add group column
result_long <- result_long %>%
  mutate(Sample_group = factor(sample_groups[Sample]))

########################################################

p <- ggplot(result_long, aes(x = Sample, y = Group)) +
  geom_point(aes(size = Category, color = Category)) +
  scale_size_manual(values = c("<10" = 2, "10-20" = 4, "20-30" = 6, "30+" = 8)) +
  scale_color_manual(values = cb_palette) +
  facet_wrap(~Sample_group, scales = "free_x", space = "free_x" , nrow = 1)+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.minor.x = element_line(color = "grey90")
  ) +
  labs(
    x = NULL,  # remove x-axis title
    y = "Group",
    size = "Category",
    color = "Category"  )


p

ggsave("Fig3_bubble_plot.pdf", plot = p, width = 8, height = 4, useDingbats = FALSE, limitsize = FALSE)
dev.off()
