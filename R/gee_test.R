#' Generalized Estimating Equations under technical replication with robust and non-robust Wald and Rao (score) tests
#'
#' @param use_geeasy When TRUE, uses `geeasy` for gee estimation, when FALSE uses `geepack`
#' @param ... Arguments that you would pass to a regular `geepack::geeglm` call. Note that for now we only have functionality for Poisson tests with log link
#'
#' @importFrom sandwich sandwich
#' @importFrom stats coef glm pnorm
#' @importFrom rlang call_match caller_env
#' @importFrom geepack geeglm
#' @import geeasy
#'
#' @examples
#' # TODO
#'
#'
#' @export
gee_test <- function(use_geeasy = TRUE, ...) {

  cl_orig <- match.call()
  cl_orig <- call_modify(cl_orig, use_geeasy = zap())

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


  ### fit the glm, ignoring warnings about non integer inputs, as we take an estimating equations mindset
  withCallingHandlers({
    gee_result <- eval(cl, envir = rlang::caller_env()) 
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "non-integer x"))
      invokeRestart("muffleWarning")
  })


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


  #### run robust score tests
  pp <- nrow(output)
  robust_score_p <- rep(NA, length = pp)

  for (p_marginal in 1:pp) {
    if (use_geeasy) {
      id <- gee_result$id
    } else {
      id <- gee_result$geese$id
    }
    robust_score_p[p_marginal] <- robust_score_test(glm_object = gee_result,
                                                    call_to_model = cl,
                                                    param = p_marginal,
                                                    id = id)
  }

  output <- cbind(output, "Robust Score p" = robust_score_p)

  output <- output[, c("Estimate","Robust Std Error", "Robust Wald p", "Robust Score p")]
  
  if (use_geeasy) {
    output <- rbind(output,
                    "correlation:alpha" = c(gee_result$geese$alpha, # "Estimate"
                                            NA, # no estimate of "Robust Std Error" from geeasy
                                            NA, # "Robust Wald p"
                                            NA)) # "Robust Score p"
  } else {
    output <- rbind(output,
                    "correlation:alpha" = c(gee_result$geese$alpha, # "Estimate"
                                            gee_result$geese$valpha.ajs, # "Robust Std Error"
                                            NA, # "Robust Wald p"
                                            NA)) # "Robust Score p"
  }
  output
}
