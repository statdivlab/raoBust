#' Generalized Linear Models with robust and non-robust Wald and Rao (score) tests
#'
#' @param ... Arguments that you would pass to a regular `glm` call. Note that for now we only have functionality for Poisson tests with log link
#'
#' @importFrom sandwich sandwich
#' @importFrom stats coef glm pnorm
#'
#' @examples
#' # glm_test(dist ~ speed, data = cars, family=poisson(link="log"))
#'
#'
#' @export
#'
glm_test <- function(...) {

  #### estimate parameters with `glm`
  glm_result <- glm(...)

  if (glm_result$family$family != "poisson") {
    stop("Amy has only implemented this for Poisson tests. Please open a GitHub issue if you're interested in other families.")
  }
  if (glm_result$family$link != "log") {
    stop("Amy has only implemented this for Poisson tests with log link. Please open a GitHub issue if you're interested in other link functions.")
  }

  output <- coef(summary(glm_result))[, -3]                 ## "z value"
  colnames(output)[3] <- "Non-robust Wald p"                ## "Pr(>|z|)"
  colnames(output)[2] <- "Non-robust Std Error"             ## "Std. Error"


  #### run robust score tests
  pp <- nrow(output)
  robust_score_ts <- vector("numeric", length = pp)
  for (p_marginal in 1:pp) {
    robust_score_ts[p_marginal] <- robust_score_test(my_formula = glm_result$formula,
                                                     df = glm_result$data,
                                                     param = p_marginal)
  }

  output <- cbind(output, "Robust Score p" = robust_score_ts)

  robust_wald_ses <- sqrt(diag(sandwich::sandwich(glm(...),
                                                  adjust=TRUE)))
  robust_wald_ps <- 2*(1-pnorm(abs(output[, "Estimate"] / robust_wald_ses)))

  output <- cbind(output, "Robust Std Error" = robust_wald_ses)
  output <- cbind(output, "Robust Wald p" = robust_wald_ps)

  output[, c("Estimate", "Non-robust Std Error", "Robust Std Error", "Non-robust Wald p", "Robust Wald p", "Robust Score p")]

}
