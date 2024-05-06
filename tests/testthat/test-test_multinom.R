nsim <- 1000
nn <- 500
jj <- 5

simulate_null_data_mult <- function(global = FALSE, ms = 10000, j = NULL) {
  beta0s <- rnorm(n = jj - 1, mean = 0, sd = 0.1)

  if (global) {
    beta1s <- rnorm(n = jj - 1, mean = 0, sd = 0)
  } else {
    beta1s <- rnorm(n = jj - 1, mean = 0, sd = 0.1) ## TODO check reasonable; make variable
  }

  if (!is.null(j)) {
    beta1s[j] <- 0
  }

  covariate1 <- cbind(1, runif(nn, 0, 1))
  xbetas <- covariate1 %*% rbind(beta0s, beta1s)
  exp_xbetas <- cbind(exp(xbetas), 1)
  probs <- exp_xbetas / rowSums(exp_xbetas)
  yy <- apply(X=probs, MARGIN=1, FUN=rmultinom, n = 1, size=ms) %>% t
  df <- list(X = matrix(covariate1[,2], ncol = 1),
             Y = yy)
  return(df)
}

df <- simulate_null_data_mult(global=TRUE)
# get_multinom_score(df$X, df$Y, joint = TRUE)


test_that("global multinomial test controls T1E", {

  ### comment out for now

  #
  # set.seed(240506)
  # ts <- rep(NA, nsim)
  # for (i in 1:nsim) {
  #   df <- simulate_null_data_mult(global=TRUE)
  #   ts[i] <- get_multinom_score(df$X, df$Y, joint = TRUE)
  # }
  # ts[ts < 0] <- 0 ### TODO investigate; throw warning
  #
  # ps <- pchisq(ts, df = jj - 1, lower.tail=FALSE)
  # hist(ps)
  #
  # qqplot(ps, runif(1000))
  #
  # 100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10)) ## should be 1, 5, 10% to be nominal
  # # for n = 500, 1000 sims: (3.6  6.4 10.5)%. Too high. Investigate.

})


df <- simulate_null_data_mult(global=FALSE, j = 2)
# get_multinom_score(df$X, df$Y, joint=FALSE, j=2)
## TODO write marginal test
