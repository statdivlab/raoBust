#' Compute score contribution to robust score statistic for glm robust score tests.
#'
#' @param i Index of observation.
#' @param model_fits The fitted glm under the null hypothesis.
#' @param yy Vector of observed responses.
#' @param xx Design matrix for model.
#' @param family The model family for the fitted glm.
#' @param link The link function utilized in the fitted glm.
#' 

score_contribution <- function(i, model_fits, yy, xx, family, link) { ### output is p x 1
  score_mat <- matrix(NA, nrow = length(model_fits), ncol = 1)
  
  if (family == "poisson" & link == "log") {
    score_vec <- matrix(model_fits[i]^2 * (yy[i] - model_fits[i]) * xx[i, ], ncol = 1)
  }
  
  if (family == "binomial" & link == "logit") {
    score_vec <- matrix(model_fits[i]*(1 - model_fits[i])^(-1)*(yy[i] - model_fits[i]) * xx[i,], ncol = 1)
  }
  
  # if (family == "gaussian" & link == "identity") {
  #   score_vec <- 
  # }
  
  return(score_vec)
}