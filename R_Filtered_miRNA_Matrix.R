# =======================================================================
#         R_miRNA_Matrix_Normalization.R
#
#   Filtering and normalization of miRNA counts matrix using edgeR 
#   package from BioConductor  
#
#   Sergio Grande García
# =======================================================================

# ===========================================================
# 1. Data setup ---------------------------------------------
# ===========================================================

# Check BioConductor repository installation 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("edgeR", quietly = TRUE))
  BiocManager::install("edgeR")
# Libraries
library(limma)
library(edgeR)
library(openxlsx)
library(plyr)
library(tidyverse)
library(ggplot2)

# ===========================================================
# 2. miRNA counts matrix preprocess -------------------------
# ===========================================================
# Database load
counts <- read.xlsx("../processed_data/allHIV_wout_outliers_norm_miRNAmatrix_cpm.txt")
raw_clinical <-  read.xlsx("../raw_data/datos_basal/BBDD_clinica_basal_preprocessed_wout_outliers.xlsx")

# Filtering of those miRNAs whose expression is higher than 10 copies
# Unnecessary step - filterBy() function of edgeR does the same process
# counts50 <- counts [which(counts$read_count>10),]

# microRNA filtering with double counts
# In the miRNA array there are repeated microRNAs, as they may come from different precursors.
# The edgeR package needs ids of each unique micro-unique. There are two options:
# i) Make the identifier unique by adding the precursor name to the id.
# (This may bias the statistical analysis because of the way miRDeep2 aligns the sequences).
# ii) Filter and select the miRNA with the most counts.


# Filtering using the method that select the miRNA with the most counts
counts_mature <- ddply(counts, .(`#miRNA`), function(x) x[which.max(x$read_count),])
row.names(counts_mature) <- counts_mature$`#miRNA`      
input_counts <- counts_mature[,5:ncol(counts_mature)]

# Subset study patiens: HIV basal moment (HIV, HIV+/HCV- Y HIV+/HCV+)
input_counts <- input_counts %>% select(starts_with(c("H", "L", "V")))

# ===========================================================
# 3. Grouping patients --------------------------------------
# ===========================================================
# We define the group to which patient belongs
group <- colnames(input_counts)
group[grepl("^V",ignore.case = T, group)] <- "VIH"
group[grepl("^L",ignore.case = T, group)] <- "VIH+/VHC-"
group[grepl("^H",ignore.case = T, group)] <- "VIH+/VHC+"

study_groups <- data.frame(group)
study_groups$group <- as.factor(study_groups$group)
rownames(study_groups) <- colnames(input_counts)
colnames(study_groups) <- "Grupos"
mycolors <- c("#F0E442", "#0072B2", "#009E73", "#D55E00")
save(study_groups,
     mycolors,
     file = "results/Study_groups_correspondence.R")


# ===========================================================
# 4. Matrix Filtering and Normalization ---------------------
# ===========================================================
# edgeR object
y <- DGEList(counts = input_counts)
y$samples[,1] <- group %>% as.factor()

# Removal of miRNAs with less than 10 reads and calculation of normalization factors (TMM)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = "TMM")
counts_df <-  as.data.frame(y$counts)
# CPM (Counts Per Million) & log scaling
counts_cpm <- cpm(y, log = F, normalized.lib.sizes = T)
pseudocounts_cpm <- cpm(y, log = T, normalized.lib.sizes = T)

# ===========================================================
# 5. Output -------------------------------------------------
# ===========================================================
#CPM
write.table(counts_cpm, 
            file = "../processed_data/allHIV_norm_miRNAmatrix_cpm.txt",
            sep = '\t', row.names = T, col.names = T, quote = T) 
#Log
write.table(pseudocounts_cpm, 
            file = "../processed_data/allHIV_norm_miRNAmatrix_pseudocounts.txt",
            sep = '\t', row.names = T, col.names = T, quote = T) 


# ===========================================================
# 6. Normalization impact visualization ---------------------
# ===========================================================
# Normalized data
pseudocounts_cpm <- cpm(y, log = T,normalized.lib.sizes = T)
boxplot(pseudocounts_cpm, xlab="Pacientes de estudio", ylab="LogCPM", 
        main = "LogCPM normalizados", xaxt = "n", las = 1)
abline(h=median(pseudocounts_cpm),col="blue")
# Not normalized data
pseudocounts_cpm_nonorm <- cpm(y, log = T,normalized.lib.sizes = F)
boxplot(pseudocounts_cpm_nonorm, xlab="Pacientes de estudio", ylab="LogCPM", 
        main = "LogCPM no normalizados", xaxt = "n", las = 1)
abline(h=median(pseudocounts_cpm_nonorm),col="blue")


