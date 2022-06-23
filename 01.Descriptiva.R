#           R_Clinical_Descriptive_All_Groups.R
#   
#   Estadística descriptiva de datos clínicos y epidemiológicos. Toma de muestra: momento = basal.
#
#         1. Todos los grupos: HC, VIH, VIH/VHC- y VIH/VHC+ 
#         2. Control sano vs pacientes VIH
#         3. Compareción pareada entre todos los grupos
#
#  SGG & DVM
#  =========================================================================
#                 01.Clinical_Descriptive.R
#  =========================================================================
#
# Descriptive statistics of clincal and epidemiological data in the basal 
# moment. Three analysis were made
#     1. All groups: stratified by HCV exposure + Healty controls (HC)
#     2. HC vs All HIV
#     3. All groups compared two by two
# For the format output Daniel Valle-Millares script tableone_nasformat was used
#
# Sergio Grande Garcia
#  =========================================================================

# =====================================================
# 0. Data setup----------------------------------------
# =====================================================
# Libraries
library(openxlsx)
library(tableone)
library(purrr)
library(tidyverse)
source("R/tableone_nasformat.R") #Daniel Valle-Millares script

# Datasets
bd <- read.xlsx("processed_data/BBDD_clinica_basal_preprocessed.xlsx") 
bd$Grupo <- bd$Grupo %>% as.factor()
bd <- bd %>% 
  dplyr::filter(., !is.na(CODIGO.CRG.basal))
bd <- bd %>% 
  dplyr::filter(., Grupo != "Control sano") %>% 
  droplevels()
load ("R/work_vars")


# ==================================================
# 1) All groups ------------------------------------
# ==================================================

# Descriptive analysis of the clinical features of patients
# stratified by HCV exposure (HIV vs. HIV+/HCV- vs. HIV+/HCV+ vs. HC)
# Descriptive table with the CreateTableOne package and modification of the 
# format with the script tableone_nasformat.R (Daniel Valle-Millares)
#
#   Sergio Grande

# Descriptive statistics w/ 'tableone'
tab1 <- CreateTableOne(vars = work_var,
                       data = bd,
                       strata = "Grupo",
                       factorVars = work_catvar)
#Table w/ numerical variables
Cont_Table <- print(tab1$ContTable,
                    showAllLevels = T,
                    nonnormal = work_nvar,
                    exact = work_catvar,
                    quote = F, noSpace = T, printToggle = T,
                    formatOptions = list(bigmark = ",", scientific = FALSE))
#Table w/ categorical variables
Cat_Table <- print(tab1$CatTable,
                   showAllLevels = T,
                   nonnormal = work_nvar,
                   exact = work_catvar,
                   quote = F, noSpace = T, printToggle = T,
                   formatOptions = list(bigmark = ",", scientific = FALSE))

# Reformat the output of categorical variables using 'nasformat()' to obtain 
# the frequency of occurrences over the total number of valid cases (!= NA)
cc <- 1
out_list <- list()
for(varname in work_catvar){
  tmp_df <- lapply(seq_along(tab1$CatTable), tab_one = tab1$CatTable,  
                   var = varname, FUN = tableone_nasformat) %>% 
    reduce(inner_join, by = "level") %>%
    mutate(variable=c(varname,rep("-", length(level)-1)), .) %>%
    relocate(variable)
  tmp_df[1,"variable"]<-varname
  out_list[[cc]]<- tmp_df
  cc <- cc + 1 
}
out_tab <- do.call(rbind, out_list)
head(out_tab)

# Export every output into an excell sheet
Export_list <- list(out_tab, Cont_Table, Cat_Table)
names(Export_list[]) <- c("Cat_nasformat_Table","Cont_Table", "Cat_Table")
#write.xlsx(Export_list, file = "results/descriptive_statistics/Descriptive_Table_All_Groups.xlsx", 
#          colNames = T, rowNames = T)

# =======================================
# 2) HC vs HIV patients (HIV + SC + CHC) 
# =======================================

# Descriptive analysis of the clinical characteristics of patients
# without differentiating by HCV exposure (HIV vs HC group)
# Descriptive table with the CreateTableOne package and modification of the 
# format with the script tableone_nasformat.R (Daniel Valle-Millares)
#
# Sergio Grande García

# New variable to compact HCV exposure groups
bd_2 <- bd %>% mutate(VIH_HC=case_when(
  Grupo == "Control sano" ~ "HC",
  Grupo != "Control sano" ~ "VIH")
  )

# Descriptive statistics w/ 'tableone'
tab2 <- CreateTableOne(vars = work_var,
                       data = bd_2,
                       strata = "VIH_HC",
                       factorVars = work_catvar)
#Table w/ numerical variables
Cont_Table_2 <- print(tab2$ContTable,
                    showAllLevels = T,
                    nonnormal = work_nvar,
                    exact = work_catvar,
                    quote = F, noSpace = T, printToggle = T,
                    formatOptions = list(bigmark = ",", scientific = FALSE))
#Table w/ categorical variables
Cat_Table_2 <- print(tab2$CatTable,
                   showAllLevels = T,
                   nonnormal = work_nvar,
                   exact = work_catvar,
                   quote = F, noSpace = T, printToggle = T,
                   formatOptions = list(bigmark = ",", scientific = FALSE))

# Reformat the output of categorical variables using 'nasformat()' to obtain 
# the frequency of occurrences over the total number of valid cases (!= NA)
cc <- 1
out_list2 <- list()
for(varname in work_catvar){
  tmp_df <- lapply(seq_along(tab2$CatTable), tab_one = tab2$CatTable,  
                   var = varname, FUN = tableone_nasformat) %>% 
    reduce(inner_join, by = "level") %>%
    mutate(variable=c(varname,rep("-", length(level)-1)), .) %>%
    relocate(variable)
  tmp_df[1,"variable"]<-varname
  out_list2[[cc]]<- tmp_df
  cc <- cc + 1 
}
out_tab2 <- do.call(rbind, out_list2)
head(out_tab2)

# Export every output into an excell sheet
Export_list <- list(out_tab2, Cont_Table_2, Cat_Table_2)
names(Export_list[]) <- c("Cat_nasformat_Table","Cont_Table", "Cat_Table")
#write.xlsx(Export_list, file = "results/descriptive_statistics/Descriptive_Table_VIH_HC.xlsx", colNames = T, rowNames = T)

# =======================================
# 3) Two by two group comparation
# =======================================

# Descriptive analysis of the clinical characteristics of patients
# comparing all study groups 2 by 2 
# Descriptive table with the CreateTableOne package and modification of the 
# format with the script tableone_nasformat.R (Daniel Valle-Millares)
#
# Sergio Grande

#================================================================
# Function definition  -----------------------------------------
#================================================================

# Function to create the 3 tables with the descriptive analysis and the inference
#Database = Work dataset 
#group = stratifying variable
#variables = All variables to include in the table
#catvar = categorical variables 
#nvar = numeric variables 
Combn_tab <- function (database, group, variables, catvar, nvar){
  tab1 <- CreateTableOne(vars = variables,
                         data = database,
                         strata = group,
                         factorVars = catvar,
                         testNonNormal = wilcox.test)
  
  res_list <- list() #List where every df is going to be added
  #Numeric variables
  res_list[[1]] <- print(tab1$ContTable,
                         showAllLevels = T,
                         nonnormal = nvar,
                         exact = catvar,
                         quote = F, noSpace = T, printToggle = T,
                         formatOptions = list(bigmark = ",", scientific = FALSE)
  )
  #Categorical variables
  res_list[[2]] <- print(tab1$CatTable,
                         showAllLevels = T,
                         nonnormal = nvar,
                         exact = catvar,
                         quote = F, noSpace = T, printToggle = T,
                         formatOptions = list(bigmark = ",", scientific = FALSE)
  )
  names(res_list) <- c("cont_tab", "cat_tab")
  return(res_list)
}

# END Functions --------------------------------
Combination <- combn(levels(bd$Grupo), 2, simplify = T) # Matrix with every permutation between study groups

# Loop to create a list containing each possible database arising from each 
# permutation of the study groups.
cc <- 1                                                              #Internal counter
list_Two_by_Two <- list()                                            #List where df is goint to be added
for (i in 1:ncol(Combination)) {
  df_temp <- subset(bd, Grupo %in% Combination[,i])                  #Subsetting the database
  df_temp$Grupo <- factor(df_temp$Grupo)                             #Restart group factor
  list_Two_by_Two[[i]] <- df_temp
  cc <- cc + 1
}                                                     

# Apply the fuction to create the descritive tables to every above databases
list_descritive_tables <- lapply(list_Two_by_Two, 
                                 function (x) Combn_tab(x, "Grupo", work_var, work_catvar, work_nvar)) 

# Excell output
# Define sheets names
xlsx.sheet.names <- c(#"Cont_table_VIH_HC", "Cat_table_VIH_HC",
                      #"Cont_table_HC_VIH+VHC-","Cat_table_HC_VIH+VHC-",
                      #"Cont_table_HC_VIH+VHC+","Cat_table_HC_VIH+VHC+",
                      "Cont_table_VIH_VIH+VHC-", "Cat_table_VIH_VIH+VHC-",
                      "Cont_table_VIH_VIH+VHC+", "Cat_table_VIH_VIH+VHC+",
                      "Cont_table_VIH+VHC+_VIH+VHC-", "Cat_table_VIH+VHC+_VIH+VHC-")

# Process table to export
# WARNING -> Do not run flatten function twice
list_descritive_tables <- purrr::flatten(list_descritive_tables)     #Match all levels in the list
names(list_descritive_tables) <- xlsx.sheet.names                    #Rename with the sheet name
for (i in 1:length(list_descritive_tables)) {
  list_descritive_tables[[i]] <- as.data.frame(list_descritive_tables[[i]])
}                                                    #Convert the objects in the list to df

# Output
write.xlsx(list_descritive_tables, file = "results/descriptive_statistics/Descriptive_Table_Two_by_Two_Groups.xlsx", 
           colNames = T, rowNames = T)

