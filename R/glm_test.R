#' Generalized Linear Models with robust and non-robust Wald and Rao (score) tests
#'
#' @param ... Arguments that you would pass to a regular `glm` call. Note that for now we only have functionality for Poisson tests with log link
#'
#' @importFrom sandwich sandwich
#' @importFrom stats coef glm pnorm
#' @import rlang
#'
#' @examples
#' # glm_test(dist ~ speed, data = cars, family=poisson(link="log"))
#'
#'
#' @export
#'
glm_test <- function(...) {

  # https://stackoverflow.com/questions/20680399/how-to-wrap-glm-in-a-function-pass-dotdotdot-directly-to-another-function-fail
  cl <- match.call()
  cl[1] <- call("glm")

  ### fit the glm, ignoring warnings about non integer inputs, as we take an estimating equations mindset
  withCallingHandlers({
    glm_result <- eval(cl, envir = rlang::caller_env()) ### is the issue here??
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "non-integer x"))
      invokeRestart("muffleWarning")
  })


  if (glm_result$family$family != "poisson") {
    stop(paste("Amy has only implemented this for Poisson families.\n",
               "You requested", glm_result$family$family, "\n",
               "Please open a GitHub issue if you're interested in other families."))
  }
  if (glm_result$family$link != "log") {
    stop(paste("Amy has only implemented this for Poisson families with log link.\n",
               "You requested link", glm_result$family$link, "\n",
               "Please open a GitHub issue if you're interested in other link functions"))
  }
  
  output <- coef(summary(glm_result))[, -3]                 ## "z value"
  colnames(output)[3] <- "Non-robust Wald p"                ## "Pr(>|z|)"
  colnames(output)[2] <- "Non-robust Std Error"             ## "Std. Error"


  # stop("no")
  #### run robust score tests
  pp <- nrow(output)
  robust_score_p <- vector("numeric", length = pp)
  for (p_marginal in 1:pp) {
    robust_score_p[p_marginal] <- robust_score_test(glm_object = glm_result,
                                                     call_to_model = cl,
                                                     param = p_marginal)
  }

  output <- cbind(output, "Robust Score p" = robust_score_p)

  robust_wald_ses <- sqrt(diag(sandwich::sandwich(glm_result,
                                                  adjust=TRUE)))
  robust_wald_ps <- 2*(1-pnorm(abs(output[, "Estimate"] / robust_wald_ses)))
  
  ci_lower <- output[,'Estimate'] - qnorm(0.975)*robust_wald_ses
  ci_upper <- output[,'Estimate'] + qnorm(0.975)*robust_wald_ses
  
  output <- cbind(output, "Robust Std Error" = robust_wald_ses)
  output <- cbind(output, "Robust Wald p" = robust_wald_ps)
  output <- cbind(output, "Lower 95% CI" = ci_lower)
  output <- cbind(output, "Upper 95% CI" = ci_upper)
  
  output <- output[, c("Estimate", "Non-robust Std Error", "Robust Std Error", "Lower 95% CI", "Upper 95% CI", "Non-robust Wald p", "Robust Wald p", "Robust Score p")]
  
  result <- list("call" = cl,
                 "coef_tab" = output)
  
  return(structure(result, class = "raoFit"))

}
