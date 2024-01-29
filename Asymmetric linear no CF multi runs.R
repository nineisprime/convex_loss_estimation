setwd("~/Desktop/Min's project/R-code/score matching new")
source("scoreFun.R")
library(MASS)

N_sim <- 1
N <- 100 # sample size
d <- 2 # dimension of X

# tail of the noise
light_tail = FALSE
############ Method A: KDE based ###############
kernel <- "gaussian"
truncation_lim <- NULL
set_to_zero <- FALSE
n_iter_kde = 40 #root finding, Newton method 
############ Method B: spline based ###############
lambda0grid = 0.1^(0:4)
robust = FALSE
n_iter_spline = 40

squared_losses <- matrix(0, nrow = N_sim, ncol = 3)

for (j in 1:N_sim){
  beta_0 <- rnorm(d)
  beta_0 <- beta_0 / sum(beta_0^2)^0.5 * 3
  mu <- 2
  
  ## Generate the data
  X <- rnorm(N * d, 0 + 1, 1)
  X <- matrix(X ,nrow = N, ncol = d)
  
  # Additive noise
  # Noise <- runif(N, -1, 1)
  # Noise <- rnorm(N, 0, 1)
  Noise <- rcauchy(N)
  
  # Asymmetric log-concave noise with mean = median = 0
  # Noise <- rnorm(N, 0, 1)
  # negNoise <- (Noise < 0)
  # Noise <- runif(N, -2*sqrt(2 / pi), 0) * negNoise + Noise * (1 - negNoise)
  
  # Noise <- rnorm(N, 0, 1)
  # negNoise <- (Noise < 0)
  # Noise1 <- rnorm(N, 0, 1)
  # Noise1 <- -abs(Noise1)^0.3
  # Noise <- Noise1 * negNoise + Noise * (1 - negNoise)
  
  # Gaussian mixture
  # N1 <- rbinom(1, N, 0.4)
  # Noise <- rnorm(N1, 3, .5)
  # Noise <- c(rnorm(N-N1, 0, 1), Noise)
  
  Y <- X %*% beta_0 + mu + Noise
  X_intercept <- cbind(rep(1, N), X)
  
  if(light_tail){
    ###use OLS as pilot
    betapilot <- solve(t(X_intercept) %*% X_intercept, t(X_intercept) %*% Y) # betapilot[1] is the intercept
  }else{
    ###use median regressor as pilot
    beta_init <- c(mu, beta_0)
    result <- optim(par = beta_init, fn = mae_loss, X = X_intercept, Y = Y, method = "L-BFGS-B")
    betapilot <- result$par
  }
  
  squared_losses[j, 1] <- sum((betapilot[-1]-beta_0)^2)
  
  Residualsorted <- sort(Y - X_intercept %*% betapilot)
  deltaknots <- diff(Residualsorted)
  
  ############# Method A: KDE based #############
  fcn_psi <- kde_decr_score_est(Residualsorted, k=1000, kernel=kernel,
                                kernel_pts=2^20, truncation_lim=truncation_lim,
                                set_to_zero=set_to_zero)
  fcn_psi_deriv <- fcn_psi[["psi_deriv"]]
  fcn_psi <- fcn_psi[["psi"]]

  betahat <- betapilot
  residuals <- Y - X_intercept %*% betahat
  psibeta <- t(X_intercept) %*% fcn_psi(residuals)
  psiprimebeta <- fcn_psi_deriv(residuals)
  psiprimebeta <- -t(X_intercept) %*% diag(psiprimebeta) %*% X_intercept
  
  for (l in 1:n_iter_kde) {
    alpha <- 1
    if(sum(psibeta^2)^0.5 < (d / N)^0.5 * 1e-9) {
      break
    }
    Update <- solve(psiprimebeta + diag(d+1), psibeta)
    # line search
    betatemp <- betahat - alpha * Update
    residuals <- Y - X_intercept %*% betatemp
    psibetatemp <- t(X_intercept) %*% fcn_psi(residuals)
    while (sum(psibetatemp^2) > sum(psibeta^2)){
      # print("shrink the step size")
      alpha <- 0.8 * alpha
      betatemp <- betahat - alpha * Update
      residuals <- Y - X_intercept %*% betatemp
      psibetatemp <- t(X_intercept) %*% fcn_psi(residuals)
    }
    betahat <- betatemp
    psiprimebeta <- fcn_psi_deriv(residuals)
    psiprimebeta <- -t(X_intercept) %*% diag(psiprimebeta) %*% X_intercept
    psibeta <- psibetatemp
  }
  if (sum(psibeta^2)^0.5 > (d + 1) * 1e-6 / N^0.5){
    print("Z equation wasn't solved properly.")
    cat("Final norm of the gradient of the loss", sum(psibeta^2)^0.5, "\n")
    squared_losses[j, 2] <- NA
  }else{squared_losses[j, 2] <- sum((betahat[-1]-beta_0)^2)}

  # ##### for inference ##### 
  # info <- mean((fcn_psi(residuals))^2)
  
  ############# Method B: Spline based, CV for lambda tuning#############
  n <- floor(N/3)
  X1 <- X[1:n, ]
  Y1 <- Y[1:n]
  X2 <- X[(n + 1):(2 * n), ]
  Y2 <- Y[(n + 1):(2 * n)]
  X3 <- X[(2 * n + 1):N, ]
  Y3 <- Y[(2 * n + 1):N]
  
  # X1_intercept <- cbind(rep(1, n), X1)
  # if(light_tail){
  #   ########## OLS as pilot  ########## 
  #   betapilot1 <- solve(t(X1_intercept) %*% X1_intercept, t(X1_intercept) %*% Y1)
  # }else{
  #   ########## median as pilot  ########## 
  #   result <- optim(par = beta_init, fn = mae_loss, X = X1_intercept, Y = Y1, method = "L-BFGS-B")
  #   betapilot1 <- result$par
  # }
  # Residualsorted1 <- sort(Y1 - X1_intercept %*% betapilot1)
  
  ########## use the global pilot  ##########
  Residualsorted1 <- sort(Y1 - X1 %*% betapilot[-1]-betapilot[1])
  cat("CV starts \n")
  MSE <- rep(0, length(lambda0grid))
  for (i in 1:length(lambda0grid)) {
    psistars <- antitonicScoreMatch(Residualsorted1, lambda0grid[i])
    psistarprimes <- psistars[["psistarprimes"]]
    psistars <- psistars[["psistars"]]
    ####reduce the memory use, not necessary
    non_redundent_index <- diff(psistars) != 0
    temp <- c(non_redundent_index, FALSE)
    non_redundent_index <- c(TRUE, non_redundent_index)
    non_redundent_index <- non_redundent_index | temp
    
    betahat2 <-spline_linear_regression_newton(beta_init, X2, Y2, Residualsorted1[non_redundent_index],
                                               psistars[non_redundent_index], 
                                               psistarprimes[non_redundent_index],
                                               n_iter = n_iter_spline, robust = robust)
    if(light_tail){
      #MSE loss for lighter tailed noises
      MSE[i] <- mean((Y3- X3 %*% betahat2[-1]-betapilot[1])^2)
    }else{
      #mean absolute error for heavy tailed noises
      X3_intercept <- cbind(rep(1, length(Y3)), X3)
      MSE[i] <- mae_loss(c(betapilot[1],betahat2[-1]), X3_intercept, Y3)
    }
  }
  lambda0 <- lambda0grid[which.min(MSE)]
  cat("CV done, lambda = ", lambda0, "\n")
  # print(MSE)
  psistars <- antitonicScoreMatch(Residualsorted, lambda0)
  psistarprimes <- psistars[["psistarprimes"]]
  psistars <- psistars[["psistars"]]

  betahat <-spline_linear_regression_newton(betapilot, X, Y, Residualsorted, psistars, psistarprimes,
                                            n_iter = n_iter_spline, robust = robust)
  
  squared_losses[j, 3] <- sum((betahat[-1]-beta_0)^2)
  
}
  
cat("Pilot MSE: ", mean(squared_losses[, 1]), '\n')
cat("kde SCALE MSE: ", mean(squared_losses[, 2]), '\n')
cat("Spline SCALE MSE: ", mean(squared_losses[, 3]), '\n')
