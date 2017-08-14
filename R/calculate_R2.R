#' Extract R^2 from a linear model
#' 
#' @description Calculate the adjusted r-squared from a linear model given
#' a dataframe and formula
#' 
#' @param data The data frame to perform the calculation on
#' @param formula A formula for the model
#' 
#' @return A scalar numeric representing the adjusted r-squared for the model
#' 
#' @author Sam Levin
#' 
#' @importFrom stats as.formula lm
#' @export
#' 

r2_calc <- function(data, formula){
  
  mod <- summary(lm(as.formula(formula), data = data))
  
  out <- mod$adj.r.squared
  
  return(out)
}
