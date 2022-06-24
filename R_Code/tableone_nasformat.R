#' @title Modify CreateTableOne output
#' @description tableone_nasformat.R edits a table of descriptive statistic (categorical variables only) obtained with tableone::CreateTableOne()$CatTable or tableone::CreateCatTable(), in such a way that non-missing values are indicated as a part of the output. 
#' @usage tableone_nasformat(
#' i,
#' tab_one,
#' var
#' )
#' @param i  Id of levels in the stratifying(grouping) variable. It can be an integer or string.
#' @param tab_one output from tableone::CatTable() or tableone::CreateTableOne()$CatTab.
#' @param var A single variable name.
#' @details No details so far.
#' @importFrom dplyr %>% lapply mutate rename
#' @import purrr reduce
#' @value A data.frame.
#' @examples 
#' # Load 
#' library(tableone)
#' library(dplyr)
#' library(purrr)
#' source(path_to/tableone_nasformat.R)
#' 
#' # Build a table of descriptive statistics with talbeone
#' tab1 <- CreateTableOne(vars = allvars,
#'                        data = bd,
#'                        strata = "Grupo",
#'                        factorVars = catvars)
#'                        
#' # Run each variable at a time                  
#' x <- lapply(seq_along(tab1$CatTable), tab_one = tab1$CatTable,  var = "SEXO", FUN = tableone_nasformat) %>%
#' reduce(inner_join, by = "level")
#' 
#' # OR execute them all
#' cc <- 1
#' out_list <- list()
#' for(varname in catvars){
#'   tmp_df <- lapply(seq_along(tab1$CatTable), tab_one = tab1$CatTable,  var = varname, FUN = tableone_nasformat) %>% 
#'     reduce(inner_join, by = "level") %>%
#'     mutate(variable=c(varname,rep("-", length(level)-1)), .) %>%
#'     relocate(variable)
#'   tmp_df[1,"variable"]<-varname
#'   out_list[[cc]]<- tmp_df
#'   cc <- cc + 1 
#' }
#' out_tab <- do.call(rbind, out_list)
#' @export
tableone_nasformat <- function(i, tab_one, var){
  gname <- names(tab_one)[i]
  dat <- tab_one[[i]][[var]]
  dat <- dat %>% 
    mutate(
      formated = paste(freq ,"/", (n-miss)," ","(",ifelse(!(percent %in% c(0,NaN)),format(round(percent, digits = 2), nsmall = 2), "0.00"),"%",")", sep = "")) %>%
    select(level,formated) %>%
    rename(!!gname := formated)
  return(dat)
}
