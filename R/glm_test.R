#' Generalized Linear Models with robust and non-robust Wald and Rao (score) tests
#'
#' @param ... Arguments that you would pass to a regular `glm` call. Any observations with `NA` values in the data (response or covariates) will be dropped.
#'
#' @importFrom sandwich sandwich
#' @importFrom stats coef glm pnorm qnorm
#' @import rlang
#'
#' @examples
#' glm_test(dist ~ speed, data = cars, family=poisson(link="log"))
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

  if (!(is.null(glm_result$na.action))) {
    message(paste0(length(glm_result$na.action), 
                   " of your observations contain missing values. These observations will be dropped from the analysis."))
  }
  
  glm_family <- glm_result$family$family
  glm_link <- glm_result$family$link
  if ((glm_family != "poisson" | glm_link != "log") & (glm_family != "binomial" | glm_link != "logit")) {
    stop(paste("This is only implemented this for Poisson family with log link and Binomial family with logit link.\n",
               "You requested", glm_family, "family and", glm_link, "link \n",
               "Please open a GitHub issue if you're interested in other families."))
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
  
  #compute robust score test for null of all coefs (except intercept) being 0
  null_model_p <- robust_score_test(glm_object = glm_result,
                                 call_to_model = cl,
                                 param = 2:pp)
  
  result <- list("call" = cl,
                 "coef_tab" = output,
                 "pval" = null_model_p)
  
  return(structure(result, class = "raoFit"))

}
