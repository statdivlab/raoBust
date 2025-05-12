#' Generalized Estimating Equations under technical replication with robust and non-robust Wald and Rao (score) tests
#'
#' @param use_geeasy When TRUE, uses `geeasy` for gee estimation, when FALSE uses `geepack`
#' @param use_jack_se When TRUE uses jackknife standard errors (which take longer), when FALSE uses sandwich standard errors
#' @param cluster_corr_coef Optional within-cluster correlation coefficient. This will only be used when parameter estimation with a GEE fails and estimation must 
#' instead be performed with a GLM. 
#' @param skip_gee When TRUE doesn't try to optimize with a GEE (just uses a GLM). This should only be used internally for testing.
#' @param ... Arguments that you would pass to a regular `geepack::geeglm` call. Note that for now we only have functionality for Poisson tests with log link
#'
#' @importFrom sandwich sandwich vcovJK
#' @importFrom stats coef glm pnorm
#' @importFrom rlang call_match caller_env
#' @importFrom geepack geeglm
#' @importFrom stats glm
#' @import geeasy
#'
#' @examples
#' # TODO
#'
#'
#' @export
gee_test <- function(use_geeasy = TRUE, use_jack_se = FALSE, cluster_corr_coef = NULL, skip_gee = FALSE, ...) {
  
  cl_orig <- match.call()
  cl_orig <- call_modify(cl_orig, use_geeasy = zap(), use_jack_se = zap(), 
                         cluster_corr_coef = zap(), skip_gee = zap())

  cl <- cl_orig
  if (use_geeasy) {
    cl[1] <- call("geelm")
  } else {
    cl[1] <- call("geeglm")
  }
  
  if (is.null(cl$id)) {
    stop("Missing replicates information (`id`). If no replicates, use `glm_test()`.")
  }

  if (!is.null(cl$weights)) {
    wts <- eval(cl$data, envir = rlang::caller_env())[[cl$weights]]
    if(!all(wts[1] == wts)) {
      warning("Amy isn't quite sure how weights interact with the replication structure. Perhaps chat with her before using?")
    }
  }

  id_col <- cl$id

  ## enforce robust Wald using "std.err" = "san.se"
  ## Documentation for geepack suggests bias can be substantial for sandwich se's with small number of clusters
  ## so choose approximate jackknife as std.err="jack" when using geepack 
  cl <- call_modify(cl,
                    "id" = as.factor(eval(cl$data, envir = rlang::caller_env())[[cl$id]]),
                    "corstr" = "exchangeable")
  if (!use_geeasy) {
    cl <- call_modify(cl, "std.err" = "jack")
  } 
  the_reorder <- order(eval(cl$data, envir = rlang::caller_env())[[id_col]])
  cl <- call_modify(cl,
                    "data" = eval(cl$data, envir = rlang::caller_env())[the_reorder, ],
                    "id" = eval(cl$id, envir = rlang::caller_env())[the_reorder])


  ### fit the gee, ignoring warnings about non integer inputs, as we take an estimating equations mindset
  withCallingHandlers({
    gee_result <- try({eval(cl, envir = rlang::caller_env())})
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "non-integer x"))
      invokeRestart("muffleWarning")
  })
  
  gee_fails <- FALSE
  if (skip_gee == TRUE) {
    gee_fails <- TRUE
  }
  if (gee_fails | inherits(gee_result, "try-error")) {
    warning("GEE solver failed. Estimation will be done with a GLM, which will provide consistent parameter estimates but will not be as efficient as estimates from a GEE. Robust standard errors will be estimated using a clustered jackknife procedure.")
    gee_fails <- TRUE
    new_cl <- cl
    new_cl[1] <- call("glm")
    new_cl <- call_modify(new_cl, corstr = zap(), id = zap())
    withCallingHandlers({
      glm_result <- try({eval(new_cl, envir = rlang::caller_env())})
    }, warning=function(w) {
      if (startsWith(conditionMessage(w), "non-integer x"))
        invokeRestart("muffleWarning")
    })
    if (inherits(glm_result, "try-error")) {
      stop("Estimation cannot be done using a GLM solver. This model cannot be optimized.")
    }
    gee_result <- glm_result
  }


  if (gee_result$family$family != "poisson") {
    stop(paste("Amy has only implemented this for Poisson families.\n",
               "You requested", gee_result$family$family, "\n",
               "Please open a GitHub issue if you're interested in other families."))
  }
  if (gee_result$family$link != "log") {
    stop(paste("Amy has only implemented this for Poisson families with log link.\n",
               "You requested link", gee_result$family$link, "\n",
               "Please open a GitHub issue if you're interested in other link functions"))
  }

  output <- coef(summary(gee_result))[, -3]                 ## "Wald"
  colnames(output)[3] <- "Robust Wald p"                ## "Pr(>|W|)"
  colnames(output)[2] <- "Robust Std Error"             ## "Std. Error"

  # get jackknife standard errors if needed (if using geeasy engine or if using glm solver)
  if (use_jack_se & use_geeasy) {
    gee_result_new <- gee_result
    gee_result_new$call$formula <- eval(gee_result$call$formula, rlang::caller_env())
    output[, 2] <- jackknife_se(object = gee_result_new, 
                                dat = eval(cl$data, envir = rlang::caller_env())[the_reorder, ],
                                id = eval(cl$id, envir = rlang::caller_env())[the_reorder])
    robust_wald_stats <- (output[, 1]/output[, 2])^2
    for (r in 1:nrow(output)) {
      output[r, 3] <- 1 - pchisq(robust_wald_stats[r], df = 1)
    }
  }
  if (gee_fails) {
    output[, 2] <- sqrt(diag(sandwich::vcovJK(glm_result, 
                                    cluster = eval(cl$id, envir = rlang::caller_env())[the_reorder])))
    robust_wald_stats <- (output[, 1]/output[, 2])^2
    for (r in 1:nrow(output)) {
      output[r, 3] <- 1 - pchisq(robust_wald_stats[r], df = 1)
    }
    
    # use user-input cluster correlation coefficient if given 
    gee_result$geese$alpha <- NULL
    if (is.null(cluster_corr_coef)) {
      warn("Because estimation has been done with a GLM, there is no estimated cluster correlation coefficient. Therefore the robust score test cannot be run. In order to run a robust score test for this model, please input a within-cluster correlation coefficient to use.")
    } else {
      gee_result$geese$alpha <- cluster_corr_coef
    }
  }
  
  ### compute 95% confidence intervals using robust std errors
  ci_lower <- output[,'Estimate'] - qnorm(0.975)*output[,'Robust Std Error']
  ci_upper <- output[,'Estimate'] + qnorm(0.975)*output[,'Robust Std Error']
  
  output <- cbind(output, "Lower 95% CI" = ci_lower)
  output <- cbind(output, "Upper 95% CI" = ci_upper)
  
  #### run robust score tests
  pp <- nrow(output)
  robust_score_p <- rep(NA, length = pp)

  if (!(gee_fails) | (gee_fails & !(is.null(cluster_corr_coef)))) {
    for (p_marginal in 1:pp) {
      if (gee_fails) {
        id <- eval(cl$id, envir = rlang::caller_env())[the_reorder]
      } else {
        if (use_geeasy) {
          id <- gee_result$id
        } else {
          id <- gee_result$geese$id
        }
      }
      robust_score_p[p_marginal] <- robust_score_test(glm_object = gee_result,
                                                      call_to_model = cl,
                                                      param = p_marginal,
                                                      id = id)
    }
  }
  
  output <- cbind(output, "Robust Score p" = robust_score_p)

  output <- output[, c("Estimate","Robust Std Error", "Lower 95% CI", "Upper 95% CI", "Robust Wald p", "Robust Score p")]
  
  if (gee_fails) {
    if (is.null(cluster_corr_coef)) {
      output <- rbind(output,
                      "correlation:alpha" = c(NA, # no estimate of correlation coefficient from GLM fit
                                              NA, # no estimate of "Robust Std Error" 
                                              NA, # no estimate of "Robust Std Error", so no CI computed
                                              NA, # no estimate of "Robust Std Error", so no CI computed
                                              NA, # "Robust Wald p"
                                              NA)) # "Robust Score p"
    } else {
      output <- rbind(output,
                      "correlation:alpha" = c(cluster_corr_coef, # user-provided cluster correlation coefficient
                                              NA, # no estimate of "Robust Std Error" 
                                              NA, # no estimate of "Robust Std Error", so no CI computed
                                              NA, # no estimate of "Robust Std Error", so no CI computed
                                              NA, # "Robust Wald p"
                                              NA)) # "Robust Score p"
    }
  } else {
    if (use_geeasy) {
      output <- rbind(output,
                      "correlation:alpha" = c(gee_result$geese$alpha, # "Estimate"
                                              NA, # no estimate of "Robust Std Error" from geeasy
                                              NA, # no estimate of "Robust Std Error", so no CI computed
                                              NA, # no estimate of "Robust Std Error", so no CI computed
                                              NA, # "Robust Wald p"
                                              NA)) # "Robust Score p"
    } else {
      output <- rbind(output,
                      "correlation:alpha" = c(gee_result$geese$alpha, # "Estimate"
                                              gee_result$geese$valpha.ajs, # "Robust Std Error",
                                              gee_result$geese$alpha - qnorm(0.95)*gee_result$geese$valpha.ajs, #95% ci lower
                                              gee_result$geese$alpha + qnorm(0.95)*gee_result$geese$valpha.ajs, #95% ci upper
                                              NA, # "Robust Wald p"
                                              NA)) # "Robust Score p"
    }
  }
  
  result <- list("call" = cl_orig,
                 "coef_tab" = output)
  return(structure(result, class = "raoFit"))
}
