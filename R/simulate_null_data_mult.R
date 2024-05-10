#' Simulate multinomial data under the null (strong or weak)
#'
#' By default, \eqn{p = 2} with a single Uniform(0,1) continuous covariate.
#'
#' @param nn Number of observations
#' @param strong Simulate under the strong null? Defaults to FALSE (simulate under the weak null)
#' @param jj Number of taxa
#' @param ms Number of counts per sample
#' @param jj_null For the weak null, which taxon should be null
#' @param sd_beta0s The beta0's are drawn from a normal distribution with mean zero. This is the standard deviation of that distribution.
#' @param sd_beta1s The beta1's are drawn from a normal distribution with mean zero. This is the standard deviation of that distribution.
#' @param overdispersion An additional normal random variable can be added to the link function to add dispersion above a multinomial distribution. This is the standard aviation for that normal variable. Useful for confirming error rate control under model misspecification.
#'
#'
#' @importFrom stats rnorm rmultinom
#'
#' @author Amy Willis
#'
#' @export
#'
simulate_null_data_mult <- function(nn,
                                    strong = FALSE,
                                    jj = 5,
                                    ms = 10000,
                                    jj_null = NULL,
                                    sd_beta0s = NULL,
                                    sd_beta1s = NULL,
                                    overdispersion = 0) {

  if (is.null(sd_beta0s)) {
    sd_beta0s <- 1
  }

  if (is.null(sd_beta1s)) {
    sd_beta1s <- 1
  }

  if (strong) {
    sd_beta1s <- 0
  }

  beta0s <- rnorm(n = jj - 1, mean = 0, sd = sd_beta0s)
  beta1s <- rnorm(n = jj - 1, mean = 0, sd = sd_beta1s)

  if (!is.null(jj_null)) {
    beta1s[jj_null] <- 0
  }

  covariate1 <- cbind(1, seq(from = 0, to = 1, length.out = nn))
  xbetas <- covariate1 %*% rbind(beta0s, beta1s)

  if (overdispersion > 0) {
    xbetas <- xbetas + rnorm(length(xbetas), mean = 0, sd = overdispersion)
  }

  exp_xbetas <- cbind(exp(xbetas), 1)
  probs <- exp_xbetas / rowSums(exp_xbetas)
  yy <- apply(X=probs, MARGIN=1, FUN=rmultinom, n = 1, size=ms) %>% t
  df <- list(X = matrix(covariate1[,2], ncol = 1),
             Y = yy)
  return(df)
}
