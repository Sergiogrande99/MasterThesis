#  ============================================================================
#                           Reservoir_by_sex_diferences.R
#  ============================================================================
#
# In this scrip a comparison of the size of the HIV reservoir is made
# in each group of patients according to sex 
# 
# Sergio Grande García
#  ============================================================================

# ====================================
# 0. Data setup ----------------------
# ====================================
#Libraries
library(tidyverse)
library(openxlsx)
library(MASS)
# Load datasets 
raw_clinical <- read.xlsx("../raw_data/datos_basal/BBDD_clinica_basal_preprocessed_wout_outliers.xlsx")
# Define every variable as their corresponding class
raw_clinical$Grupo <-  raw_clinical$Grupo %>% as.factor()
raw_clinical$SEXO <-  raw_clinical$SEXO %>% as.factor()
raw_clinical$reservorio.CD4r.basal <- raw_clinical$reservorio.CD4r.basal %>% as.numeric()

# =================================================================
# Data summary of HIV reservoir - Strata by HCV exposure and sex
# =================================================================
#Simple way to visualize in All HIV (No HCV exposure strata)
summary(raw_clinical$reservorio.CD4r.basal[raw_clinical$SEXO == "HOMBRE"])
summary(raw_clinical$reservorio.CD4r.basal[raw_clinical$SEXO == "MUJER"])
summary(raw_clinical$reservorio.total[raw_clinical$SEXO == "HOMBRE"])
summary(raw_clinical$reservorio.total[raw_clinical$SEXO == "MUJER"])

# Console printing of every data summary output - Strata by HCV exposure and sex
# Matrix with every permutation between HCV exposure group and sex
Combination <- expand.grid(levels(raw_clinical$Grupo), levels(raw_clinical$SEXO)) %>% as.matrix()
for (i in 1:6) { #Loop to define output and know each permutation output
 df.temp <-  dplyr::filter(raw_clinical, raw_clinical$Grupo %in% Combination[i,1] & raw_clinical$SEXO %in% Combination[i,2])
   df.temp$reservorio.CD4r.basal <- df.temp$reservorio.CD4r.basal %>% as.numeric()
   print(paste("Resumen de la variable Reservorio CD4r en pacientes", Combination[i,1], "y", Combination[i,2], sep = " "))
   print(summary(df.temp$reservorio.CD4r.basal))
   print(paste("Resumen de la variable Reservorio total en pacientes", Combination[i,1], "y", Combination[i,2], sep = " "))
   print(summary(df.temp$reservorio.total))
}

# =================================================================
# Differences in CD4r HIV reservoir by HCV exposure and sex
# =================================================================

wilcox.test(raw_clinical$reservorio.CD4r.basal ~ raw_clinical$SEXO)
# Monoinfected
clinical <- raw_clinical[raw_clinical$Grupo == "VIH+",]
wilcox.test(clinical$reservorio.CD4r.basal ~ clinical$SEXO)
# Spontaneous clearance
clinical <- raw_clinical[raw_clinical$Grupo == "VIH+VHC-",]
wilcox.test(clinical$reservorio.CD4r.basal ~ clinical$SEXO)
# Coinfected
clinical <- raw_clinical[raw_clinical$Grupo == "VIH+VHC+",]
wilcox.test(clinical$reservorio.CD4r.basal ~ clinical$SEXO)

# =================================================================
# Differences in total PBMCs HIV reservoir by HCV exposure and sex
# =================================================================

wilcox.test(raw_clinical$reservorio.total ~ raw_clinical$SEXO)
# Monoinfected
clinical <- raw_clinical[raw_clinical$Grupo == "VIH+",]
wilcox.test(clinical$reservorio.total ~ clinical$SEXO)
# Spontaneous clearance
clinical <- raw_clinical[raw_clinical$Grupo == "VIH+VHC-",]
wilcox.test(clinical$reservorio.total ~ clinical$SEXO)
# Coinfected
clinical <- raw_clinical[raw_clinical$Grupo == "VIH+VHC+",]
wilcox.test(clinical$reservorio.total ~ clinical$SEXO)







