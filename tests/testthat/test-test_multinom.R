nsim <- 200
nn <- 100
jj <- 5
simulate_null_data_mult <- function(global = FALSE) {
  beta0s <- rnorm(n = jj - 1, mean = 0, sd = 0.1)

  if (global) {
    beta1s <- rnorm(n = jj - 1, mean = 0, sd = 0)
  } else {
    beta1s <- rnorm(n = jj - 1, mean = 0, sd = 0.1)
  }
  covariate1 <- cbind(1, runif(nn, 0, 1))
  xbetas <- covariate1 %*% rbind(beta0s, beta1s)
  exp_xbetas <- cbind(exp(xbetas), 1)
  probs <- exp_xbetas / rowSums(exp_xbetas)
  yy <- apply(X=probs, MARGIN=1, FUN=rmultinom, n = 1, size=10000) %>% t
  df <- list(X = matrix(covariate1[,2], ncol = 1),
             Y = yy)
  return(df)
}

df <- simulate_null_data_mult()
get_multinom_score(df$X, df$Y, joint = FALSE, j = 2)


test_that("global multinomial test controls T1E", {

  ts <- rep(NA, nsim)
  for (i in 1:nsim) {
    df <- simulate_null_data_mult(global=TRUE)
    ts[i] <- get_multinom_score(df$X, df$Y, joint = TRUE)
  }
  ts[ts < 0] <- 0
  summary(ts)

  pchisq(12, df = jj - 1, lower.tail=FALSE)

  qqplot(ts, rchisq(1000, 4))

  ps <- pchisq(ts, df = jj - 1, lower.tail=FALSE)
  hist(ps)
  100*c(mean(ps < 0.01), mean(ps < 0.05), mean(ps < 0.10)) ## should be 1, 5, 10% to be nominal

})
