# Wald confidence interval of robustlmm::rlmer beta coefficients

# Adapted from code provided Ben Bolker on StackExchange: https://stats.stackexchange.com/questions/233800/how-can-i-get-confidence-intervals-for-fixed-effects-using-the-rlmer-function-r

confint.rlmerMod <- function(object, level = 0.95) {
  
  # Extract beta coefficients
  beta <- fixef(object)
  
  # Extract names of coefficients
  parm <- names(beta)
  
  # Extract standard errors for the coefficients
  se <- sqrt(diag(vcov(object)))
  
  # Set level of confidence interval
  z <- qnorm((1 + level) / 2)
  
  # Calculate CI
  ctab <- cbind(beta - (z * se), 
                beta + (z * se))
  
  # label column names
  colnames(ctab) <- c(paste(100 * ((1 - level) / 2), '%'),
                      paste(100 * ((1 + level) / 2), '%'))
  
  # Output
  return(ctab[parm, ])
}