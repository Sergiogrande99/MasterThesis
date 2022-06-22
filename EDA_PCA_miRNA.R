# ===================================================================
#       EDA_PCA_miRNA.R
# ===================================================================
#
# Exploratory data analysis to find clusters and see how patients 
# (all HIV+) are grouped according to miRNA expression profile.
# We will look for possible differences according to:
# 1. Sex 
# 2. HCV exposure
# 
# HCV exposure will be stratified by HCV exposure and groupings will 
# be sought according to:
# 1. Sex
#
# Sergio Grande Garcia
# ===================================================================

# ===================================================================
# 0) Data setup -----------------------------------------------------
# ===================================================================

#Libraries
library(ggplot2)
library(scales)
library(grid)
library(ggbiplot)
library(pca3d)
library(tidyverse)
library(openxlsx)
library(factoextra)
library(rgl) # gif library
#Load data 
cpm_counts <- read.table("../processed_data/allHIV_norm_miRNAmatrix_cpm.txt",
                         sep = '\t', header = T)
cpm_counts <- t(cpm_counts)
load("R/Study_groups_correspondence.R")
study_groups <- study_groups[rownames(study_groups) != "HVN" & 
                               rownames(study_groups) != "LVE" & 
                               study_groups$Grupos != "Control sano", ,drop = F]
clinical <- read.xlsx("../processed_data/BBDD_clinica_basal_preprocessed_wout_outliers.xlsx")
# Subset patiens (Only patients on both datasets)
rownames(clinical) <- clinical$CODIGO.CRG.basal
clinical <- clinical %>%  
  dplyr::filter(., CODIGO.CRG.basal %in% rownames(cpm_counts)) 

 
# ============================================================
# 1) Analisis superficial exploratorio -----------
# ============================================================

# PCA object
cpm_counts_pca <- prcomp(cpm_counts, center = TRUE, scale. = TRUE) 
#First visualization
plot(cpm_counts_pca$x[,1], cpm_counts_pca$x[,2])             
#Scree plot - variance explanation (Eigenvalue)
scree_mirna <- fviz_eig(cpm_counts_pca, addlabels = TRUE, ylim = c(0, 50), 
         xlab =  "Dimensiones",
         ylab = "Porcentaje de contribución a la varianza", 
         main = "Contribución de cada dimensión a la dispersión de los datos",
         ggtheme = theme_classic())   

#Otros
fviz_pca_var(cpm_counts_pca, col.var = "black")                   #Eigenvector of each vector
fviz_cos2(cpm_counts_pca, choice = "var", axes = 1:2)             #Representation quality
fviz_contrib(cpm_counts_pca, choice = "var", axes = 1, top = 10)  #Contribution of each variable to PC1

##K-means clusters
counts_transform = as.data.frame(cpm_counts_pca$x)
fviz_nbclust(cpm_counts, kmeans, method = 'wss')
fviz_nbclust(cpm_counts, kmeans, method = 'silhouette')
fviz_nbclust(cpm_counts, kmeans, method = 'gap_stat')

kmeans(cpm_counts, centers = 2, nstart = 50)
fviz_cluster(., data = cpm_counts)


# ============================================================
# 2) Graphs -------------------------------------------------
# ============================================================

#2D
# By HCV Exposure
ggbiplot(cpm_counts_pca, var.axes = F, 
         groups = clinical$Grupo,        # Delete groups argument to visualize PCA without grouping
         choices = c(1,2)) +             # Choices argument changes dimensions of the analysis
  ggtitle("Análisis de Componentes Principales",
          "Diferenciación de clusters en función del perfil de expresión de miRNA") +
  theme_classic()+
  xlab(paste("Componente Principal 1 (", format(round(scree_mirna$data$eig[1], 2), nsmall = 2), " % Varianza)"))+
  ylab(paste("Componente Principal 2 (", format(round(scree_mirna$data$eig[2], 2), nsmall = 2), " % Varianza)"))+
  labs(color="Grupos") +
  scale_color_manual(name="Grupos", values= mycolors[-1])
# By Sex
ggbiplot(cpm_counts_pca, var.axes = F, groups = clinical$SEXO, 
         ellipse = T, choices = c(1,2)) +             #Choices argument changes dimensions of the analysis
  ggtitle("Análisis de Componentes Principales",
          "Diferenciación de clusters en función del perfil de expresión de miRNA") +
  theme_classic()+
  xlab(paste("Componente Principal 1 (", format(round(scree_mirna$data$eig[1], 2), nsmall = 2), " % Varianza)"))+
  ylab(paste("Componente Principal 2 (", format(round(scree_mirna$data$eig[2], 2), nsmall = 2), " % Varianza)"))+
  labs(color="Sexo") +
  scale_color_manual(name="Sexo", values= c("blue", "red"))
  
##3D
pca3d(cpm_counts_pca, group = study_groups$Grupos,   # Delete groups argument to visualize PCA without grouping
      show.ellipses=TRUE, show.plane=FALSE, legend = "topleft", 
      components = c(1,3,2), #Components argument changes dimensions of the analysis
      palette = mycolors) 
dir.create("animation_merge")

#GIF
for (i in 1:360) {
  view3d(userMatrix=rotationMatrix(2*pi * i/360, 0, 1, 0))
  rgl.snapshot(filename=paste("animation_merge/frame-",
                              sprintf("%03d", i), ".png", sep=""))}

#Source: https://2-bitbio.com/2017/04/animated-3d-pca-plots-in-r.html



                        