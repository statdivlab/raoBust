#' Generalized Estimating Equations under technical replication with robust and non-robust Wald and Rao (score) tests
#'
#' @param ... Arguments that you would pass to a regular `geepack::geeglm` call. Note that for now we only have functionality for Poisson tests with log link
#'
#' @importFrom sandwich sandwich
#' @importFrom stats coef glm pnorm
#' @importFrom rlang call_match caller_env
#' @importFrom geepack geeglm
#'
#' @examples
#' TODO
#'
#'
#' @export
gee_test <- function(...) {

  # cl_orig <- rlang::call_match(defaults = TRUE)
  cl_orig <- match.call()

  cl <- cl_orig
  cl[1] <- call("geeglm")

  if (is.null(cl$id)) {
    stop("Missing replicates information (`id`). If no replicates, use `glm_test()`.")
  }

  if (!is.null(cl$weights)) {
    wts <- eval(cl$data)[[cl$weights]]
    if(!all(wts[1] == wts)) {
      warning("Amy isn't quite sure how weights interact with the replication structure. Perhaps chat with her before using?")
    }
  }

  # print(cl)

  id_col <- cl$id


  # print(eval(cl$data))
  # print(cl$id)
  # print(eval(cl$data)[[cl$id]])


  ## enforce robust Wald using "std.err" = "san.se"
  ## Documentation for geepack suggests bias can be substantial for sandwich se's with small # of clusters
  ## so choose approximate jackknife as std.err="jack"
  cl <- call_modify(cl,
                    "id" = as.factor(eval(cl$data)[[cl$id]]),
                    "corstr" = "exchangeable",
                    "std.err" = "jack")
  the_reorder <- order(eval(cl$data)[[id_col]])
  cl <- call_modify(cl,
                    "data" = eval(cl$data)[the_reorder, ],
                    "id" = eval(cl$id)[the_reorder])


  ### fit the glm, ignoring warnings about non integer inputs, as we take an estimating equations mindset
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

    # robust_score_p[p_marginal] <- robust_score_test(glm_object = geeglm_result,
    #                                                 call_to_model = cl,
    #                                                 param = p_marginal)
    robust_score_p[p_marginal] <- NA
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
