# ====================================================================
#          AnalisisExplo_datosClinicos.R
#   
# Exploratory descriptive analysis of HIV reservoir. In this script 
# we explored the variables in order to find outliers, describe the 
# variable distribution and visualize the data
#
# Sergio Grande Garc�a & Daniel Valle-Millares
# ====================================================================

#Libraries
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(EnvStats)
library(ggforce)
library(ggdist)
library(gghalves)
#====================================================================
# 0. Data setup -----------------------------------------------------
#====================================================================
bd <- read.xlsx("processed_data/BBDD_clinica_basal_preprocessed.xlsx")
bd <- bd %>% mutate(VIH_HC=case_when( Grupo == "Control sano" ~ "HC", Grupo != "Control sano" ~ "VIH"))
bd$reservorio.CD4r.basal <- as.numeric(bd$reservorio.CD4r.basal)
bd_infected <- bd %>% subset(Grupo != "Control sano")

# ==========================
# 1. Reservoir vs. sex -----
# ==========================
# ggplot2 visualization
# HIV reservoir in CD4r
cd4r_sexo <- ggplot(bd_infected, aes(y = reservorio.CD4r.basal, x = Grupo, fill = SEXO)) +
  geom_boxplot(outlier.alpha = 0.5) +
  theme(legend.position="bottom", plot.title = element_text(face="bold"), text = element_text(size = 20)) + 
  labs(title="Reservorio CD4r",  fill = "Sexo", 
       x ="Exposición al VHC", y = "Reservorio VIH (copias/millón cel.)")+
  scale_fill_discrete(name = "Sexo", labels = c("Hombre", "Mujer"))

# HIV reservoir in PBMCs CD4-
restos_sexo <- ggplot(bd_infected, aes(y =  reservorio.Restos.basal, x = Grupo, fill = SEXO)) +
  geom_boxplot(outlier.alpha = 0.5) + 
  theme(legend.position="bottom",plot.title = element_text(face="bold"), text = element_text(size = 20))+
  labs(title="Reservorio CMSPs CD4-",  fill = "Sexo",
       x ="Exposición al VHC", y = "Reservorio VIH (copias/millón cel.)")+
  scale_fill_discrete(name = "Sexo", labels = c("Hombre", "Mujer"))

# HIV reservoir in total PBMCs
restotal_sexo <- ggplot(bd_infected, aes(y =  reservorio.total, x = Grupo, fill = SEXO)) +
  geom_boxplot(outlier.alpha = 0.5) + 
  theme(legend.position="bottom", plot.title = element_text(face="bold"), text = element_text(size = 20))+
  labs(title="Reservorio total",  fill = "Sexo",
       x ="Exposición al VHC", y = "Reservorio VIH (copias/millón cel.)")+
  scale_fill_discrete(name = "Sexo", labels = c("Hombre", "Mujer"))


#Representation
plot_grid(cd4r_sexo, restos_sexo, restotal_sexo, 
          ncol = 3, nrow = 1,label_size = 10, label_y = 1.03)

# =================================================
# 2. Analisis de outliers - Rosner's test ---------
# =================================================
# Hypothesis testing
t_cd4r      <- rosnerTest(bd_infected$reservorio.CD4r.basal, k = 3)
t_restos    <- rosnerTest(bd_infected$reservorio.Restos.basal, k = 2)
t_restotal  <- rosnerTest(bd_infected$reservorio.total, k = 3)

# Output process
res_cd4r <- t_cd4r$all.stats %>% mutate(Variable = "Reservorio CD4r", Paciente = slice(.data = bd_infected, Obs.Num)$CODIGO.CRG.basal)
res_restos <- t_restos$all.stats %>% mutate(Variable = "Reservorio restos", Paciente = slice(.data = bd_infected, Obs.Num)$CODIGO.CRG.basal)
res_restotal <- t_restotal$all.stats %>% mutate(Variable = "Reservorio total", Paciente = slice(.data = bd_infected, Obs.Num)$CODIGO.CRG.basal)
outlier_results <- rbind(res_cd4r, res_restos, res_restotal)

# Output export to excel
# write.xlsx(outlier_results, file = "results/01.analisis_exploratorio/Outliers_RosnersTest.xlsx")


# ggplot2 visualization
bd_infected_outfilt <- bd_infected %>% subset(!(CODIGO.CRG.basal %in% outlier_results$Paciente))

# HIV reservoir in CD4r
outfilt_cd4r_sexo <- ggplot(bd_infected_outfilt, aes(y = reservorio.CD4r.basal, x = Grupo, fill = SEXO)) +
  geom_boxplot() + theme(legend.position="bottom")
# HIV reservoir in PBMCs CD4-
outfilt_restos_sexo <- ggplot(bd_infected_outfilt, aes(y =  reservorio.Restos.basal, x = Grupo, fill = SEXO)) +
  geom_boxplot() + theme(legend.position="bottom")
# HIV reservoir in total PBMCs
outfilt_restotal_sexo <- ggplot(bd_infected_outfilt, aes(y =  reservorio.total, x = Grupo, fill = SEXO)) +
  geom_boxplot() + theme(legend.position="bottom")

# Exporting the data
pdf(file = "results/01.analisis_exploratorio/figuras/AE_reservorio.pdf", width = 10, height = 12)
plot_grid( cd4r_sexo, 
           restos_sexo, 
           restotal_sexo,
           outfilt_cd4r_sexo, 
           outfilt_restos_sexo, 
           outfilt_restotal_sexo ,
           labels=c("Reservorio CD4 restos", 
                    "Reservorio restos", 
                    "Reservorio total",
                    "Reservorio CD4 restos  \n(sin outliers)", 
                    "Reservorio restos \n(sin outliers)", 
                    "Reservorio total \n(sin outliers)"), ncol = 3, nrow = 2,label_size = 10)
dev.off()



# ====================================
# Outlier elimination from databases
# ====================================
# Dataset load and outlier elimination from Clinical database
bd <- read.xlsx("../raw_data/datos_basal/BBDD_Clinica_basal.xlsx")

clinical_wout_outliers <- bd[bd$CODIGO.CRG.basal != "HVN" & bd$CODIGO.CRG.basal != "LVE",]
# Output
#write.xlsx(clinical_wout_outliers, "../processed_data/BBDD_clinica_basal_preprocessed_wout_outliers.xlsx")

# Dataset load and outlier elimination from counts matrix
cpm_counts <- read.table("../processed_data/allHIV_norm_miRNAmatrix_cpm.txt")
cpm_counts_wout_outliers <- cpm_counts[, colnames(cpm_counts) != "HVN" & colnames(cpm_counts) != "LVE"]
# Output
#write.table(cpm_counts_wout_outliers, "../processed_data/allHIV_wout_outliers_norm_miRNAmatrix_cpm.txt",
           sep = '\t', row.names = T, col.names = T, quote = T)


# ====================================
# EDA w/out outliers -----------------
# ====================================
# HIV reservoir in CD4r
cd4_raincloud <- ggplot(bd_infected_outfilt, aes(y = reservorio.CD4r.basal, x = Grupo, fill = SEXO)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) +
  geom_jitter (alpha = .3, size = 2, position=position_jitterdodge(jitter.width = 0.05,
                                                         jitter.height = 0,
                                                         dodge.width = 0.2,)) +
  labs(title="Reservorio CD4r", fill = "Sexo",
       x ="Exposición al VHC", y = "Reservorio VIH (copias/millón cel.)")+
  theme(text = element_text(size = 20), legend.position = "none")
# HIV reservoir in total PBMCs
total_raincloud <- ggplot(bd_infected_outfilt, aes(y = reservorio.total, x = Grupo, fill = SEXO)) + 
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .2, outlier.shape = NA) +
  geom_jitter (alpha = .3, size = 2, position=position_jitterdodge(jitter.width = 0.05,
                                                         jitter.height = 0,
                                                         dodge.width = 0.2,)) +
  labs(title="Reservorio total", fill = "Sexo",
       x ="Exposición al VHC", y = "Reservorio VIH (copias/millón cel.)")+
  theme(text = element_text(size = 20)) +
  scale_fill_discrete(name = "Sexo", labels = c("Hombre", "Mujer"))
# Representation
plot_grid(cd4_raincloud, total_raincloud)













