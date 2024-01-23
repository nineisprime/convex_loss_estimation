

source("scoreFun.R")
library(quantreg)

N_sim <- 1
N <- 1000 # sample size
d <- 10 # dimension of X



for (j in 1:N_sim){
  beta_0 <- rnorm(d)
  beta_0 <- beta_0 / sum(beta_0^2)^0.5 * 3
  mu <- 2

  X <- rnorm(N * d, 0 + 1, 1)
  X <- matrix(X ,nrow = N, ncol = d)

  Noise <- rcauchy(N)

  Y <- X %*% beta_0 + mu + Noise

  Xtemp <- cbind(rep(1, N), X)
  res = rq.fit(x=Xtemp, y=Y, tau=0.5)

  betapilot = res$coefficients

  Residualsorted = sort(Y - Xtemp %*% betapilot)

  truncation_lim = NULL
  set_to_zero = FALSE
  fcn_psi <- kde_decr_score_est(Residualsorted, k=1000, kernel="gaussian",
                                kernel_pts=2^21, truncation_lim=truncation_lim, set_to_zero=set_to_zero)
}


