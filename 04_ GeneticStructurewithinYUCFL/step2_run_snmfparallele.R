library(LEA)
library(tidyverse)
library(hues)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
K <- as.numeric(args[1])
project_name <- paste0("/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/02_AD1FL_n166/02_FLn196_LEA/LEA/snmf_results/snmf_K", K)
setwd(project_name)

geno_file <- "./FLn196.geno"
snmf_project=snmf(geno_file, K=K, entropy= TRUE, repetitions=10, project= "new", CPU = 20)

save(snmf_project, file = paste0("snmf_K",K."RData"))
