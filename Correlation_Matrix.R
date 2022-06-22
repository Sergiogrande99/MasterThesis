# ===================================================================
#       EDA_PCA_miRNA.R
# ===================================================================
#
# Correlation analysis of HIV reservoir and expresion of each miRNA
#
# Sergio Grande Garcia
# ===================================================================

# ===================================================================
# 0. Data setup -----------------------------------------------------
# ===================================================================
# Libraries
library(openxlsx)
library(ggcorrplot)
library(corrr)
library(Hmisc)
library(RcmdrMisc)
library(corrplot)
library(mdatools)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
# Load data
raw_cpm_counts <- read.table("processed_data/allHIV_norm_miRNAmatrix_cpm.txt",
                         header = T, sep = '\t', row.names = 1)

raw_reservoir <- read.table("processed_data/Basal_Reservoir_Preprocess.txt", 
                            sep = "\t")
load("R/reservoir_variables")

# Subset patiens (Only patients in both datasets) 
reservoir <- filter(raw_reservoir, rownames(raw_reservoir) %in% colnames(raw_cpm_counts))
cpm_counts <- raw_cpm_counts[,colnames(raw_cpm_counts) %in% rownames(raw_reservoir)]
cpm_counts <- as.data.frame(t(cpm_counts))

# dataframe for correlation matrix
input_correlation <- left_join(rownames_to_column(cpm_counts), 
                               rownames_to_column(reservoir), by = "rowname")
rownames(input_correlation) <- input_correlation$rowname
input_correlation <- input_correlation[,-1]
# ===================================================================
# 1. Pearson method -------------------------------------------------
# ===================================================================
# Correlation matrix object
pval_cor <- rcorr.adjust(input_correlation, type = "pearson", use = "complete.obs")

corR <- pval_cor$R$r[rownames(pval_cor$R$r) %in% colnames(reservoir),
                     colnames(pval_cor$R$r) %in% colnames(cpm_counts)]  # Correlation coefficients matrix
corP <- pval_cor$R$P[rownames(pval_cor$R$P) %in% colnames(reservoir),
                     colnames(pval_cor$R$P) %in% colnames(cpm_counts)]  # pvalue matrix

# Filtering by correlation coefficient (<- 0.5 & > 0.5) & pvalue (p<0.05)
sig_element <- ifelse(corP < 0.05, corR, NA)
sig_element[!(sig_element >= 0.5 | sig_element <= -0.5 )] <- NA
sig_element <- sig_element[,colSums(is.na(sig_element))<nrow(sig_element)]


# ===================================================================
# 2. Spearman method -------------------------------------------------
# ===================================================================
# If Pearson method has been compiled first, run:
# rm(list = ls())
# And code in 0.Data setup

# Correlation matrix object
pval_cor <- rcorr.adjust(input_correlation, type = "spearman", use = "complete.obs")

corR <- pval_cor$R$r[rownames(pval_cor$R$r) %in% colnames(reservoir),
                     colnames(pval_cor$R$r) %in% colnames(cpm_counts)]  # Correlation coefficients matrix
corP <- pval_cor$R$P[rownames(pval_cor$R$P) %in% colnames(reservoir),
                     colnames(pval_cor$R$P) %in% colnames(cpm_counts)]  # pvalue matrix

# Filtering by correlation coefficient (<- 0.5 & > 0.5) & pvalue (p<0.05)
sig_element <- ifelse(corP < 0.05, corR, NA)
sig_element[!(sig_element >= 0.5 | sig_element <= -0.5 )] <- NA
sig_element <- sig_element[,colSums(is.na(sig_element))<nrow(sig_element)]

# ===================================================================
# 3. Output ---------------------------------------------------
# ===================================================================
write.table(sig_element, file = "processed_data/all_HIV_correlation_matrix.txt",
            sep = "\t", row.names = T, col.names = T)

# ===================================================================
# 4. Heatmap plot ---------------------------------------------------
# ===================================================================
# Filter out miRNAS with higher variability 
var_genes_CPM <- apply(as.data.frame(t(cpm_counts)), 1, var)          # Variability estimation
select_var_tw <- names(sort(var_genes_CPM, decreasing=TRUE))[1:50]    # Subseting 50 most variable - Change value for more or less miRNAs
highly_variable_tw <- (cpm_counts)[,select_var_tw]

correlation_df_test <- left_join(rownames_to_column(highly_variable_tw), rownames_to_column(reservoir), by = "rowname")
rownames(correlation_df_test) <- correlation_df_test$rowname
correlation_df_test <- correlation_df_test[,-1]
# Create correlation matrix
tw_test <- cor(correlation_df_test, use = "complete.obs")
tw_test <- tw_test[res_vars, 1:50]
#Corplot
corrplot(tw_test[res_vars,], method = "color", 
         tl.cex = 1, cl.cex = 1, tl.col = "black",
         sig.level = 0.05)
#Heatmap
pheatmap(tw_test, cluster_cols = T, color = hcl.colors(50, "RdYlBu"), 
         scale = "none", main = "Heatmap de la matriz de correlacion", 
         show_colnames = T, cluster_rows = F, show_rownames = T, border_color = NA,
         legend_breaks = c(1, 0, -1),legend_labels = c("1", "0", "-1"), display_numbers = TRUE)


