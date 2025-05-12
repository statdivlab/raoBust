#' Print function
#'
#' @param x Object of class \code{raoFit}
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{NULL}. Displays printed model summary.
#'
#' @author Shirley Mathur 
#'
#' @method print raoFit
#'
#' @export

print.raoFit <- function (x, ...) {
  
  #print out call
  cat("\nCall:\n",
     paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  
  
  cat("\nCoefficient estimates:\n")
  print(data.frame(x$coef_tab, check.names = FALSE))
  
}
