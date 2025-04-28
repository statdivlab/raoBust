#' Robust score (Rao) tests for multinomial regression.
#'
#' @param x Object returned from \code{multinom_test}
#'
#'@return \code{NULL}. Displays printed model summary.
#'
#' @author Shirley Mathur
#'
#' @export

print_multinom <- function (x) {
  
  #print out call
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  #get dimensions of coefficients
  p <- nrow(x$mle1)
  J <- ncol(x$mle1)
  
  #make table of coefficients
  coef_tab <- data.frame("category" = rep(1:J, each = p),
                         "covariate" = rep(1:p, J),
                         "estimate" = rep(NA, p*J),
                         "se" = rep(NA, p*J))
  
  #set covariate terms appropriately if formula was used in call
  if (!is.null(x$call$formula)) {
    
    #get names of variables used for model fit, and then use these to populate covariate column of output table
    coef_names <- c(labels(terms(as.formula(x$call$formula), data = get(x$call$data))))
    coef_tab$covariate <- rep(c("(intercept)",coef_names), J)
    
  }
  
  #populate estimate and se columns of output table with mle under alternative and robust wald se, respectively
  coef_tab$estimate <- c(x$mle1)
  coef_tab$se <- c(x$wald_se)
  
  cat("\nCoefficient estimates under the alternative:\n")
  print(coef_tab)
  
}
