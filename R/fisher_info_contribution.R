#' Compute fisher information contribution to robust score statistic for glm robust score tests.
#'
#' @param i Index of observation.
#' @param model_fits The fitted glm under the null hypothesis.
#' @param yy Vector of observed responses.
#' @param xx Design matrix for model.
#' @param m Number of parameters fixed under null hypothesis.
#' @param family The model family for the fitted glm.
#' @param link The link function utilized in the fitted glm.
#' 

fisher_info_contribution <- function(i, model_fits, yy, xx, m = 1, family, link) { ### output is p x p... Di^T Vi Di
  fisher_info_mat <- matrix(NA, nrow = length(model_fits), ncol = 1)
  if (family == "poisson" & link == "log") {
    fisher_info_mat <- model_fits[i] * xx[i, ] %*% t(xx[i, ])
  }
  
  if (family == "binomial" & link == "logit") {
    fisher_info_mat <- (1-model_fits[i])^(-1)*(model_fits[i])^(-1) * xx[i, ] %*% t(xx[i, ])
  }
  
  if (family == "gaussian" & link == "identity") {
    n <- length(yy)
    p <- ncol(xx)
    mod_df <- p-m
    sigma2_tilde <- sum((yy - model_fits)^2)/(n - mod_df)
    fisher_info_mat <- (1/sigma2_tilde)*xx[i, ] %*% t(xx[i, ])
  }
  
  if (family == "gaussian" & link == "log") {
    n <- length(yy)
    p <- ncol(xx)
    mod_df <- p-m
    sigma2_tilde <- sum((yy - model_fits)^2)/(n - mod_df)
    fisher_info_mat <- ((model_fits[i]^2)/sigma2_tilde)*xx[i, ] %*% t(xx[i, ])
  }
  
  return(fisher_info_mat)
}