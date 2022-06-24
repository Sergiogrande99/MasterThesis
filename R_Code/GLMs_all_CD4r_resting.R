#  =========================================================================
#                            GLMs_all_CD4r_resting.R
#  =========================================================================
#
# Association analysis of HIV reservoir size and miRNA expression profile
#
# Sergio Grande Garcia
#  =========================================================================

# ===================================================================
# 0. Data setup -----------------------------------------------------
# ===================================================================
#Libraries
library(foreign)
library(ggplot2)
library(MASS)
library(tidyverse)
library(ggiraphExtra)
library(openxlsx)

#Load data 
raw_cpm_counts <- read.table("../processed_data/allHIV_norm_miRNAmatrix_cpm.txt",
                             header = T, sep = '\t', row.names = 1)
raw_clinical <- read.xlsx("../raw_data/datos_basal/BBDD_clinica_basal_preprocessed_wout_outliers.xlsx")
load("R/Study_groups_correspondence.R")

# Subset patientes (Only patients in both datasets)
cpm_counts <- as.data.frame(t(raw_cpm_counts))
cpm_counts <- cpm_counts %>% 
  dplyr::filter(., rownames(cpm_counts) %in% raw_clinical$CODIGO.CRG.basal)    
clinical <- raw_clinical %>% 
  dplyr::filter(., CODIGO.CRG.basal %in% rownames(cpm_counts))         
rownames(clinical) <- clinical$CODIGO.CRG.basal
study_groups <- study_groups %>%
  dplyr::filter(., rownames(study_groups) %in% rownames(cpm_counts))   


# Subset clincal variables for the analysis and the model fitting
clinical_vars <- c("SEXO", "EDAD.(anos)", "Meses.de.infeccion.por.VIH", "RIESGO", 
                   "Tto.VIH.basal_clasificacion", "linfocitos.T.CD4.(celulas/mm3).(basal)",
                   "nadir.CD4.(celulas/mm3).(basal)", "IL28b.(lab)","reservorio.CD4r.basal",               
                   "reservorio.Restos.basal", "Unspliced.CD4r.basal.MEDIA", "Single.spliced.CD4r.basal.MEDIA",   
                   "Multiple.spliced.CD4r.basal.MEDIA","Unspliced.Restos.basal.MEDIA",       
                   "Single.spliced.Restos.basal.MEDIA", "Multiple.spliced.Restos.basal.MEDIA",
                   "reservorio.total")
clinical <- clinical %>%
  dplyr::select(., one_of(clinical_vars))

# Input dataframe for the model
input_glm <- left_join(rownames_to_column(study_groups), 
                       rownames_to_column(clinical),by = "rowname")

input_glm <- input_glm %>% rename(., paciente = rowname)
input_glm$Grupos <- input_glm$Grupos %>% as.factor() %>% droplevels()

# Adaptation of the counts matrix for the model
cpm_counts <- cpm_counts %>% rownames_to_column() 
cpm_counts <- cpm_counts %>% rename(., paciente = rowname)

# ===================================================================
# 1. Generalized Linear Model testing -------------------------------
# ===================================================================
# 1.1. Previous defined function for the model ----------------------
source("R/Multiple_dependent_variables_GLM.R")

# 1.2) GLM fitting --------------------------------------------------------
# Model adjusted by HCV exposure, sex, time of HIV infection, CD4 counts and HIV tto.
GLM_all_vars_tto <- glm_multiple_depend (input_df = input_glm, mir_df = cpm_counts, 
                                         formula = reservorio.total ~ mir + Grupos + SEXO + Tto.VIH.basal_clasificacion +
                                           Meses.de.infeccion.por.VIH + `linfocitos.T.CD4.(celulas/mm3).(basal)`)
GLM_all_vars_tto_out_tab <- GLM_all_vars_tto[[1]]
GLM_all_vars_tto_sig <- dplyr::filter(GLM_all_vars_tto_out_tab, GLM_all_vars_tto_out_tab$FDR <= 0.05)
# Model adjusted by HCV exposure, sex, time of HIV infection and CD4 counts
GLM_all_vars <- glm_multiple_depend (input_df = input_glm, mir_df = cpm_counts, 
                                     formula = reservorio.total ~ mir + Grupos + SEXO + 
                                       Meses.de.infeccion.por.VIH + `linfocitos.T.CD4.(celulas/mm3).(basal)`)
GLM_all_vars_out_tab <- GLM_all_vars[[1]]
GLM_all_vars_sig <- dplyr::filter(GLM_all_vars_out_tab, GLM_all_vars_out_tab$FDR <= 0.05)

# Model adjusted by HCV exposure, sex and time of HIV infection
GLM_no_CD4 <- glm_multiple_depend (input_df = input_glm, mir_df = cpm_counts, 
                                   formula = reservorio.total ~ mir + Grupos + SEXO + Meses.de.infeccion.por.VIH)
estimate_tab <- GLM_no_CD4[[1]]
GLM_no_CD4_sig <- dplyr::filter(estimate_tab, estimate_tab$FDR <= 0.05, estimate_tab$Estimate >= 1.1 | estimate_tab$Estimate <= 0.9)

# Model adjusted by HCV exposure and sex
GLM_grupo_sexo <- glm_multiple_depend (input_df = input_glm, mir_df = cpm_counts, 
                                       formula = reservorio.total ~ mir + Grupos + SEXO)
GLM_grupo_sexo_out_tab <- GLM_grupo_sexo[[1]]
GLM_grupo_sexo_sig <- dplyr::filter(GLM_grupo_sexo_out_tab, GLM_grupo_sexo_out_tab$FDR <= 0.05)

# Model adjusted by HCV exposure
GLM_grupo <- glm_multiple_depend (input_df = input_glm, mir_df = cpm_counts, 
                                  formula = reservorio.total ~ mir + Grupos)
GLM_grupo_out_tab <- GLM_grupo[[1]]
GLM_grupo_sig <- dplyr::filter(GLM_grupo_out_tab, GLM_grupo_out_tab$FDR <= 0.05)

# Not adjusted model
GLM_no_adjust <- glm_multiple_depend (input_df = input_glm, mir_df = cpm_counts, 
                                      formula = reservorio.total ~ mir)
GLM_no_adjust_out_tab <- GLM_no_adjust[[1]]
GLM_no_adjust_sig <- dplyr::filter(GLM_no_adjust_out_tab, GLM_no_adjust_out_tab$FDR <= 0.05)


# 1.3) Summary ------------------------------------------------------------
# Dataframe containing every AIC of each model
AIC_summary <- data.frame(GLM_all = as.numeric(GLM_all_vars_out_tab$AIC),
                          GLM_all_no_CD4 = as.numeric(estimate_tab$AIC),
                          GLM_HCV_Sex = as.numeric(GLM_grupo_sexo_out_tab$AIC),
                          GLM_HCV = as.numeric(GLM_grupo_out_tab$AIC),
                          GLM_no_adjust = as.numeric(GLM_no_adjust_out_tab$AIC),
                          GLM_all_vars_tto = as.numeric(GLM_all_vars_tto_out_tab$AIC))
AIC_test <- stack(AIC_summary)
# Comparison of each fitted model (AIC mean)
summary(AIC_test)
# Hypothesis testing
pairwise.t.test(AIC_test$values, AIC_test$ind, p.adjust.method = "fdr")
# ===================================================================
# 2. Selected Generalized Linear Model-------------------------------
# ===================================================================
# To ensure the biological relevance of the results, we defined an estimate threshold 
GLM_HIV_CD4r <- glm_multiple_depend (input_df = input_glm, mir_df = cpm_counts, 
                                           formula = reservorio.CD4r.basal ~ mir + Grupos + SEXO + Meses.de.infeccion.por.VIH)
GLM_HIV_CD4r_out_tab <- GLM_HIV_CD4r[[1]]
GLM_HIV_CD4r_sig <- dplyr::filter(GLM_HIV_CD4r_out_tab, GLM_HIV_CD4r_out_tab$FDR <= 0.05, 
                                  GLM_HIV_CD4r_out_tab$Estimate >= 1.1 | GLM_HIV_CD4r_out_tab$Estimate <= 0.9) #Estimate threshold

#Sex Interaction
#Male
input_glm_Male <- input_glm %>% dplyr::filter(., SEXO == "HOMBRE") 
GLM_HIV_CD4r_M <- glm_multiple_depend (input_df = input_glm_Male, mir_df = cpm_counts, 
                                             formula = reservorio.CD4r.basal ~ mir + Grupos + Meses.de.infeccion.por.VIH)
GLM_HIV_CD4r_M_out_tab <- GLM_HIV_CD4r_M[[1]]
GLM_HIV_CD4r_M_sig <- dplyr::filter(GLM_HIV_CD4r_M_out_tab, GLM_HIV_CD4r_M_out_tab$FDR <= 0.05, 
                                    GLM_HIV_CD4r_M_out_tab$Estimate >= 1.1 | GLM_HIV_CD4r_M_out_tab$Estimate <= 0.9)

#Female
input_glm_Female <- input_glm %>% dplyr::filter(., SEXO == "MUJER")
GLM_HIV_CD4r_F <- glm_multiple_depend (input_df = input_glm_Female, mir_df = cpm_counts, 
                                             formula = reservorio.CD4r.basal ~ mir + Grupos + Meses.de.infeccion.por.VIH)
GLM_HIV_CD4r_F_out_tab <- GLM_HIV_CD4r_F[[1]]
GLM_HIV_CD4r_F_sig <- dplyr::filter(GLM_HIV_CD4r_F_out_tab, GLM_HIV_CD4r_F_out_tab$FDR <= 0.05, 
                                    GLM_HIV_CD4r_F_out_tab$Estimate >= 1.1 | GLM_HIV_CD4r_F_out_tab$Estimate <= 0.9)

# QQ-plot for quality assesment
# All HIV w/out sex stratification 
miRNA_sig <- GLM_HIV_CD4r_sig$miR
GLM_statistics <- GLM_HIV_CD4r[[2]]
GLM_statistics <- GLM_statistics[miRNA_sig]
#All HIV Male
miRNA_sig_M <- GLM_HIV_CD4r_M_sig$miR
GLM_statistics_M <- GLM_HIV_CD4r_M[[2]]
GLM_statistics_M <- GLM_statistics_M[miRNA_sig_M]
# All HIV Female
miRNA_sig_F <- GLM_HIV_CD4r_F_sig$miR
GLM_statistics_F <- GLM_HIV_CD4r_F[[2]]
GLM_statistics_F <- GLM_statistics_F[miRNA_sig_F]

# ===================================================================
# 3. Output ---------------------------------------------------------
# ===================================================================
#Excell with the GLM dataframes
# All HIV w/out sex stratification 
write.xlsx(GLM_HIV_CD4r_sig, file = "results/GLMs/GLM_HIV_all_CD4r_woutOutliers.xlsx", 
           colNames = T, rowNames = F)
#All HIV Male
write.xlsx(GLM_HIV_CD4r_M_sig, file = "results/GLMs/GLM_HIV_all_Male_CD4r_woutOutliers.xlsx", 
           colNames = T, rowNames = F)
# All HIV Female
write.xlsx(GLM_HIV_CD4r_F_sig, file = "results/GLMs/GLM_HIV_all_Female_CD4r_woutOutliers.xlsx", 
           colNames = T, rowNames = F)

##PDFs with qqnorm plots
# All HIV w/out sex stratification 
pdf('results/GLMs/Normal_qqplot_GLM_res_CD4r_residues.pdf')
pdf.options(width = 9, height = 7)
for (i in 1:length(GLM_statistics)){
  print(plot(GLM_statistics[[i]], which = 2,  main = names(GLM_statistics)[[i]], sub = ""))
}
dev.off()
#All HIV Male
pdf('results/GLMs/Normal_qqplot_GLM_res_CD4r_M_residues.pdf')
pdf.options(width = 9, height = 7)
for (i in 1:length(GLM_statistics_M)){
  print(plot(GLM_statistics_M[[i]], which = 2,  main = names(GLM_statistics_M)[[i]], sub = ""))
}
dev.off()
# All HIV Female
pdf('results/GLMs/Normal_qqplot_GLM_res_CD4r_F_residues.pdf')
pdf.options(width = 9, height = 7)
for (i in 1:length(GLM_statistics_F)){
  print(plot(GLM_statistics_F[[i]], which = 2,  main = names(GLM_statistics_F)[[i]], sub = ""))
}
dev.off()




