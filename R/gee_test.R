#' Generalized Estimating Equations under technical replication with robust and non-robust Wald and Rao (score) tests
#'
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
gee_test <- function(...) {

  cl_orig <- match.call()

  cl <- cl_orig
  cl[1] <- call("geeglm")

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
  ## so choose approximate jackknife as std.err="jack"
  cl <- call_modify(cl,
                    "id" = as.factor(eval(cl$data, envir = rlang::caller_env())[[cl$id]]),
                    "corstr" = "exchangeable",
                    "std.err" = "jack")
  the_reorder <- order(eval(cl$data, envir = rlang::caller_env())[[id_col]])
  cl <- call_modify(cl,
                    "data" = eval(cl$data, envir = rlang::caller_env())[the_reorder, ],
                    "id" = eval(cl$id, envir = rlang::caller_env())[the_reorder])


  ### fit the glm, ignoring warnings about non integer inputs, as we take an estimating equations mindset
  # first, run geeasy::geelm() to see if it converges 
  # (to avoid an infinite loop with geepack::geeglm())
  cl_check <- cl
  cl_check[1] <- call("geelm")
  cl_check <- call_modify(cl_check, std.err = zap())
  geelm_result <- try({
    eval(cl_check, envir = rlang::caller_env())
  })
  if (inherits(geelm_result, "try-error")) {
    stop("GEE does not converge. Sarah and Amy are working on figuring out why!")
  }
  
  withCallingHandlers({
    geeglm_result <- eval(cl, envir = rlang::caller_env()) ### is the issue here??
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "non-integer x"))
      invokeRestart("muffleWarning")
  })


  if (geeglm_result$family$family != "poisson") {
    stop(paste("Amy has only implemented this for Poisson families.\n",
               "You requested", geeglm_result$family$family, "\n",
               "Please open a GitHub issue if you're interested in other families."))
  }
  if (geeglm_result$family$link != "log") {
    stop(paste("Amy has only implemented this for Poisson families with log link.\n",
               "You requested link", geeglm_result$family$link, "\n",
               "Please open a GitHub issue if you're interested in other link functions"))
  }

  output <- coef(summary(geeglm_result))[, -3]                 ## "Wald"
  colnames(output)[3] <- "Robust Wald p"                ## "Pr(>|W|)"
  colnames(output)[2] <- "Robust Std Error"             ## "Std. Error"


  #### run robust score tests
  pp <- nrow(output)
  robust_score_p <- rep(NA, length = pp)

  for (p_marginal in 1:pp) {

    robust_score_p[p_marginal] <- robust_score_test(glm_object = geeglm_result,
                                                    call_to_model = cl,
                                                    param = p_marginal,
                                                    id = geeglm_result$geese$id)
  }

  output <- cbind(output, "Robust Score p" = robust_score_p)

  output <- output[, c("Estimate","Robust Std Error", "Robust Wald p", "Robust Score p")]

  output <- rbind(output,
                  "correlation:alpha" = c(geeglm_result$geese$alpha, # "Estimate"
                                          geeglm_result$geese$valpha.ajs, # "Robust Std Error"
                                          NA, # "Robust Wald p"
                                          NA)) # "Robust Score p"
  output
}
