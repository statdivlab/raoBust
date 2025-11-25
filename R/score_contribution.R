#' Compute score contribution to robust score statistic for glm robust score tests.
#'
#' @param i Index of observation.
#' @param model_fits The fitted glm under the null hypothesis.
#' @param yy Vector of observed responses.
#' @param xx Design matrix for model.
#' @param m Number of parameters fixed under null hypothesis.
#' @param family The model family for the fitted glm.
#' @param link The link function utilized in the fitted glm.
#' 

score_contribution <- function(i, model_fits, yy, xx, m = 1, family, link) { ### output is p x 1
  score_mat <- matrix(NA, nrow = length(model_fits), ncol = 1)
  
  if (family == "poisson" & link == "log") {
    score_vec <- matrix((yy[i] - model_fits[i]) * xx[i, ], ncol = 1)
  }
  
  if (family == "binomial" & link == "logit") {
    #score_vec <- matrix(model_fits[i]*(1 - model_fits[i])^(-1)*(yy[i] - model_fits[i]) * xx[i,], ncol = 1)
    score_vec <- matrix((yy[i] - model_fits[i]) * xx[i, ], ncol = 1)
  }
  
  if (family == "gaussian" & link == "identity") {
    n <- length(yy)
    p <- ncol(xx)
    mod_df <- p-m
    sigma2_tilde <- sum(yy - model_fits)/(n - mod_df)
    score_vec <- matrix((1/sigma2_tilde)*(yy[i] - model_fits[i])*xx[i,], ncol = 1)
  }
  
  if (family == "gaussian" & link == "log") {
    n <- length(yy)
    p <- ncol(xx)
    mod_df <- p-m
    sigma2_tilde <- sum(yy - model_fits)/(n - mod_df)
    score_vec <- matrix((1/sigma2_tilde)*model_fits[i]*(yy[i] - model_fits[i])*xx[i,], ncol = 1)
  }

  return(score_vec)
}