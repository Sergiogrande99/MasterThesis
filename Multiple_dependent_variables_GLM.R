#  ==========================================================================
#                   Multiple_dependent_variables_GLM.R
#  ==========================================================================
#
# In this scrip I have developed a function that generate, via a loop, 
# a table that summarizes the output (miRNA, estimate, interquartile 
# range, pvalue, AIC and the adjusted pvalue (FDR)) of multiples GLMs  
# changing the independent (x) variable
#
# Sergio Grande Garcia
#  ==========================================================================

# input_df = dataframe with HIV reservoir and clinical variables
# mir_df = counts miRNA matrix 
# formula = GLM formula, must be "depend_var(y)  ~ independent_var(x) + covar"

glm_multiple_depend <- function(input_df, mir_df, formula) {
  out_list <- list()   #List containing every selected output if the GLM 
  out_glm <- list()    #List containg the whole GLM output
  cc <- 1              #Internal counter
  for (i in 2:ncol(mir_df)) {
    print(colnames(mir_df)[i])                                            #Allows to monitor the progress of the multiple models
    df_temp <- left_join(input_df, mir_df[, c(1, i)], by = "paciente")    #In each cycle it introduces a miRNA into the input df
    names(df_temp)[ncol(df_temp)] <- 'mir'                                #renames the miRNA ID in order to define just one formula
    
    glm_temp <- tryCatch({                                                #Sometimes, the model can fail so to avoid that, I used tryCatch
      glm(formula, data = df_temp, family = Gamma(link = "log"))          #Model arguments. In my case the dependent variable follows a gamma distribution and used the link function "log"
    }, error = function(e) {                                            
      "divergence errors"
    })
    if(glm_temp == "divergence errors"){                     #If the function fails, tryCatch executes this code
      mirname = colnames(mir_df)[i]                          #In my case, the output must be "ERROR"
      est= "ERROR"
      low = "ERROR"
      up = "ERROR"
      pval = "ERROR"
      aic = "ERROR"
      
    }else{                                                   #If the function works, tryCatch executes this code In my case, I define the output
      mirname = colnames(mir_df)[i]                          #miRNA analised
      est= summary(glm_temp)$coefficients[2,1] %>% exp()     #The estimate of the model
      low = tryCatch({                                       #Low quantile
        confint(glm_temp)[2,1] %>% exp()                     #Sometimes, quantile coefficients my fail to be calculated, to avoid it, tryCatch
      }, error = function(e) {
        "profiling errors"
      })
      up = tryCatch({                                        #Up quantile
        confint(glm_temp)[2,2] %>% exp()
      }, error = function(e) {
        "profiling errors"
      })
      pval = summary(glm_temp)$coefficients[2,4]            #pvalue
      aic = glm_temp$aic                                    #AIC
    }
    
    df_out <- data.frame(miR = mirname,                     #Dataframe with every value defined above
                         Estimate = est,                 
                         low_IC = low,
                         up_IC = up,
                         P_value = pval,
                         AIC = aic)
    out_list[[cc]] <- df_out                                #The dataframe is added to the list in the position of the internal counter
    out_glm[[cc]] <- glm_temp                               #The model is added to the list in the position of the internal counter
    names(out_glm)[cc] <- mirname
    cc <- cc + 1                                            #Plus 1 to the counter
  }
  out_tab <- do.call(rbind,out_list)                        #Compacting every dataframe into one
  out_tab$FDR <- p.adjust(p = out_tab$P_value, method = 'fdr')    #pvalue adjustment (FDR)
  return(list(out_tab, out_glm))
}
