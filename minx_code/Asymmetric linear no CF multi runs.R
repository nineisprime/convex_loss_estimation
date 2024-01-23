source("scoreFun.R")
library(MASS)

mae_loss <- function(beta, X, Y) {
  predictions <- X %*% beta
  mean_abs_error <- mean(abs(predictions - Y))
  return(mean_abs_error)
}

N_sim <- 1
N <- 1000 # sample size
d <- 2 # dimension of X


############ Method A: KDE based ###############
kernel <- "gaussian"
truncation_lim <- NULL
set_to_zero <- FALSE
############ Method B: spline based ###############
lambda0grid = 0.1^(0:4)
robust = FALSE
n_iter_spline = 35

squared_losses <- matrix(0, nrow = N_sim, ncol = 4)

for (j in 1:N_sim){
  beta_0 <- rnorm(d)
  beta_0 <- beta_0 / sum(beta_0^2)^0.5 * 3
  mu <- 2
  
  ## Generate the data
  X <- rnorm(N * d, 0 + 1, 1)
  X <- matrix(X ,nrow = N, ncol = d)
  
  # Additive noise
  # Noise <- runif(N, -1, 1)
  Noise <- rnorm(N, 0, 1)
  # Noise <- rcauchy(N)
  
  # Asymmetric log-concave noise with mean = median = 0
  # Noise <- rnorm(N, 0, 1)
  # negNoise <- (Noise < 0)
  # Noise <- runif(N, -2*sqrt(2 / pi), 0) * negNoise + Noise * (1 - negNoise)

  # Noise1 <- rnorm(N, 0, 1)
  # Noise1 <- -abs(Noise1)^0.3
  # Noise <- Noise1 * negNoise + Noise * (1 - negNoise)
  
  # Gaussian mixture
  # N1 <- rbinom(1, N, 0.4)
  # Noise <- rnorm(N1, 3, .5)
  # Noise <- c(rnorm(N-N1, 0, 1), Noise)
  
  
  Y <- X %*% beta_0 + mu + Noise
  
  ###use OLS as pilot
  # Xtemp <- cbind(rep(1, N), X)
  # betapilot <- solve(t(Xtemp) %*% Xtemp, t(Xtemp) %*% Y) # betapilot[1] is the intercept
  
  ###use median regressor as pilot
  beta_init <- c(mu, beta_0)
  Xtemp <- cbind(rep(1, N), X)
  result <- optim(par = beta_init, fn = mae_loss, X = Xtemp, Y = Y, method = "L-BFGS-B")
  betapilot <- result$par
  
  squared_losses[j, 1] <- sum((betapilot[-1]-beta_0)^2)
  
  Residualsorted <- sort(Y - Xtemp %*% betapilot)
  
  ############# Method A #############
  n_iter = 10 # number of Newton iterations
  fcn_psi <- kde_decr_score_est(Residualsorted, k=1000, kernel="gaussian",
                                kernel_pts=2^21, truncation_lim=truncation_lim, set_to_zero=set_to_zero)
  fcn_psi_deriv <- fcn_psi[["psi_deriv"]]
  fcn_psi <- fcn_psi[["psi"]]
  m<-1e5
  xs<-(1:m)/(m*100)-0.5
  derivs<-fcn_psi_deriv(xs)
  plot(xs,derivs)
  plot(fcn_psi, xlim = c(-4,4))
  plot(fcn_psi_deriv, xlim = c(-4,4))
  plot(fcn_psi)
  plot(fcn_psi_deriv)
  betahat <- betapilot
  residuals <- Y - Xtemp %*% betahat
  psibeta <- t(Xtemp) %*% fcn_psi(residuals)
  psiprimebeta <- fcn_psi_deriv(residuals)
  psiprimebeta <- -t(Xtemp) %*% diag(psiprimebeta) %*% Xtemp

  for (l in 1:n_iter) {
    print("newton step")
    alpha <- 1
    if(sum(psibeta^2)^0.5 < (d + 1) * 1e-6 / N^0.5) {
      break
    }

    psiprime_pseudo <- ginv(psiprimebeta) # low dimension
    Update <- psiprime_pseudo %*% psibeta
    betatemp <- betahat - alpha * Update
    # line search
    residuals <- Y - Xtemp %*% betatemp
    psibetatemp <- t(Xtemp) %*% fcn_psi(residuals)
    # psiprimebetatemp <- fcn_psi_deriv(residuals)
    # psiprimebetatemp <- -t(Xtemp) %*% diag(psiprimebetatemp) %*% Xtemp

    while (sum(psibetatemp^2) > 2 * sum(psibeta^2)){
      print("shrink the step size")
      alpha <- 0.8 * alpha
      betatemp <- betahat - alpha * Update
      residuals <- Y - Xtemp %*% betatemp
      psibetatemp <- t(Xtemp) %*% fcn_psi(residuals)
      # psiprimebetatemp <- fcn_psi_deriv(residuals)
      # psiprimebetatemp <- -t(Xtemp) %*% diag(psiprimebetatemp) %*% Xtemp
    }
    betahat <- betatemp
    psiprimebeta <- fcn_psi_deriv(residuals)
    psiprimebeta <- -t(Xtemp) %*% diag(psiprimebeta) %*% Xtemp
    psibeta <- psibetatemp
  }
  if (sum(psibeta^2)^0.5 > (d + 1) * 1e-6 / N^0.5){
    cat("Final l2 norm of the Z function value ", sum(psibeta^2)^0.5, "\n")
  }
  squared_losses[j, 2] <- sum((betahat[-1]-beta_0)^2)

  # ############# Method B with CV#############
  # n <- floor(N/3)
  # X1 <- X[1:n, ]
  # Y1 <- Y[1:n]
  # X2 <- X[(n + 1):(2 * n), ]
  # Y2 <- Y[(n + 1):(2 * n)]
  # X3 <- X[(2 * n + 1):N, ]
  # Y3 <- Y[(2 * n + 1):N]
  # 
  # ########## OLS as pilot
  # # X1temp <- cbind(rep(1, n), X1)
  # # betapilot1 <- solve(t(X1temp) %*% X1temp, t(X1temp) %*% Y1) # betapilot1[1] is the intercept
  # # Residualsorted2 <- sort(Y2 - X2 %*% betapilot1[-1]-betapilot1[1])
  # 
  # ########## median as pilot
  # X1temp <- cbind(rep(1, n), X1)
  # result <- optim(par = beta_init, fn = mae_loss, X = X1temp, Y = Y1, method = "L-BFGS-B")
  # betapilot1 <- result$par
  # Residualsorted2 <- sort(Y2 - X2 %*% betapilot1[-1]-betapilot1[1])
  # 
  # cat("CV starts \n")
  # MSE <- rep(0, length(lambda0grid))
  # # MSE1 <- rep(0, length(lambda0grid))
  # for (i in 1:length(lambda0grid)) {
  #   psistars <- antitonicScoreMatch(Residualsorted2, lambda0grid[i])
  #   psistarprimes <- psistars[["psistarprimes"]]
  #   psistars <- psistars[["psistars"]]
  # 
  #   betahat1 <-spline_linear_regression_newton(beta_init, X1, Y1, Residualsorted2, psistars, psistarprimes,
  #                                              n_iter = n_iter_spline, robust = robust)
  # 
  #   ########## OLS
  #   # X3temp <- cbind(rep(1, length(Y3)), X3)
  #   # betapilot3 <- solve(t(X3temp) %*% X3temp, t(X3temp) %*% Y3)
  # 
  #   ########## median
  #   X3temp <- cbind(rep(1, length(Y3)), X3)
  #   result <- optim(par = beta_init, fn = mae_loss, X = X3temp, Y = Y3, method = "L-BFGS-B")
  #   betapilot3 <- result$par
  # 
  #   MSE[i] <- sum((X3 %*% (betahat1[-1]-betapilot3[-1]))^2)
  #   # betahat2 <-spline_linear_regression_newton(beta_init, X1, Y1, Residualsorted1, psistars, psistarprimes,
  #   #                                            n_iter = n_iter_spline, robust = !robust, deltaknots = deltaknots1)
  #   # MSE1[i] <- sum((X2 %*% (betahat2[-1]-betapilot2[-1]))^2)
  # }
  # lambda0 <- lambda0grid[which.min(MSE)]
  # cat("CV ends, lambda = ", lambda0, "\n")
  # print(MSE)
  # psistars <- antitonicScoreMatch(Residualsorted, lambda0)
  # psistarprimes <- psistars[["psistarprimes"]]
  # psistars <- psistars[["psistars"]]
  # 
  # ########################################test section
  # # psistars <- antitonicScoreMatch(Residualsorted, 0.001)
  # # psistarprimes <- psistars[["psistarprimes"]]
  # # psistars <- psistars[["psistars"]]
  # # plot(psistars, type = "l", col = "blue", lwd = 2, main = "spline")
  # # plot(psistarprimes, type = "l", col = "blue", lwd = 2, main = "spline")
  # ########################################
  # 
  # betahat <-spline_linear_regression_newton(betapilot, X, Y, Residualsorted, psistars, psistarprimes,
  #                                           n_iter = n_iter_spline, robust = robust)
  # 
  # # lambda01 <- lambda0grid[which.min(MSE1)]
  # # psistars1 <- antitonicScoreMatch(Residualsorted, lambda01)
  # # psistarprimes1 <- psistars1[["psistarprimes"]]
  # # psistars1 <- psistars1[["psistars"]]
  # #
  # # betahat1 <-spline_linear_regression_newton(betapilot, X, Y, Residualsorted, psistars1, psistarprimes1,
  # #                                           n_iter = n_iter_spline, robust = !robust, deltaknots = deltaknots)
  # 
  # squared_losses[j, 3] <- sum((betahat[-1]-beta_0)^2)
  # # squared_losses[j, 4] <- sum((betahat1[-1]-beta_0)^2)
  
}
  
# cat("OLS MSE: ", mean(squared_losses[, 1]), '\n')
cat("Median regressor MSE: ", mean(squared_losses[, 1]), '\n')
cat("kde SCALE MSE: ", mean(squared_losses[, 2]), '\n')
cat("Spline SCALE MSE: ", mean(squared_losses[, 3]), '\n')
cat("robust Spline SCALE MSE: ", mean(squared_losses[, 4]), '\n')
