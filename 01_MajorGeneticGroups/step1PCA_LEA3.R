library(LEA)
library(tidyverse)
library(hues)
library(parallel)
library(ggplot2)

setwd(getwd())

file1 <- list.files(pattern = "\\.ped$")

ped2geno(file1)
ped2lfmm(file1)

file2 <- list.files(pattern = "\\.lfmm$")

pca.scree.hap.plot <- pca(file2,
                          scale=T, K=40) %>%
  tracy.widom %>%
  ggplot(aes(N, percentage)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels=scales::label_percent()) +
  labs(title = "PCA Scree Plot",
       x = "Principal Components",
       y = "Percentage of Variance") +
  theme_minimal()


pdf("LEA1_pca_subset2.pdf", width = 10, height = 10)
#plotcross-entropycriterionforallrunsinthesnmfproject 
pca.scree.hap.plot
dev.off()

