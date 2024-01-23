source("scoreFun.R")

N <- 1000# sample size
d <- 2 # dimension of X
beta_0 <- 3 * runif(d) + 1
mu <- 2

## Generate the data
X <- rnorm(N * d, 0, 1)
X <- array(X, c(N, d))

# Additive noise
# Noise <- runif(N, -1, 1)
# Noise <- rnorm(N, 0, 1)
Noise <- rcauchy(N)

# Asymmetric log-concave noise with mean = median = 0
# Noise <- rnorm(N, 0, 1)
# negNoise <- (Noise < 0)
# Noise <- runif(N, -2*sqrt(2 / pi), 0) * negNoise + Noise * (1 - negNoise)

# Noise1 <- rnorm(N, 0, 1)
# Noise1 <- -abs(Noise1)^0.3
# Noise <- Noise1 * negNoise + Noise * (1 - negNoise)


Y <- X %*% beta_0 + mu + Noise

###use OLS as pilot
Xtemp <- cbind(rep(1, N), X)
betapilot <- solve(t(Xtemp) %*% Xtemp, t(Xtemp) %*% Y) # betapilot[1] is the intercept
# betapilot[1] <- betapilot[1] + 3
Residualsorted <- sort(Y - Xtemp %*% betapilot)


############ Method A ###############
kernel <- "gaussian"
truncation_lim <- NULL
set_to_zero <- FALSE
############ Method B ###############
robust = FALSE

############# Method A #############
n_iter = 8
fcn_psi <- kde_decr_score_est(Residualsorted, k=1000, kernel=kernel,
                              kernel_pts=2^20, truncation_lim=truncation_lim, set_to_zero=set_to_zero)

fcn_psi_deriv <- fcn_psi[["psi_deriv"]]
fcn_psi <- fcn_psi[["psi"]]
betahat <- betapilot
Xtemp <- cbind(rep(1, N), X)
for (l in 1:n_iter) {
  residuals <- Y - Xtemp %*% betahat
  
  psibeta <- fcn_psi(residuals)
  psibeta <- t(psibeta) %*% (Xtemp)
  psiprimebeta <- fcn_psi_deriv(residuals)
  psiprimebeta <- -t(Xtemp) %*% diag(psiprimebeta) %*% Xtemp
  if (kappa(psiprimebeta) > 1e9) {
    cat("Hessian H is singular. Using H + I instead.\n")
    psiprimebeta <- psiprimebeta + diag(nrow(psiprimebeta))
  }
  # Update <- solve(psiprimebeta, t(psibeta))
  if(sum(psibeta^2)^0.5 < (d + 1) * 1e-9) {
    break
  }
  
  Update<- tryCatch({
    solve(psiprimebeta, t(psibeta))
  }, error = function(err) {
    print("H is singular; H <- H + I")
    solve(psiprimebeta + diag(nrow(psiprimebeta)), t(psibeta))
  })
  
  # Updatenorm <- sum(Update^2)
  # cat("Update's Euclidean norm ", Updatenorm^0.5, '\n')
  # if(Updatenorm < 1e-17) {
  #   betahat <- betahat-Update
  #   break
  # }
  betahat <- betahat-Update
}
###check if betahat solve the Z-equation
residuals <- Y - Xtemp %*% betahat
psibeta <- fcn_psi(residuals)
psibeta <- t(psibeta) %*% (Xtemp)
cat("psibetahat has squared Euclidean norm ", sum(psibeta^2), '\n')
cat("KDE based SCALE's squared error ", sum((betahat[-1]-beta_0)^2), '\n')
cat("OLS's squared error ", sum((betapilot[-1]-beta_0)^2), '\n')



############# Method B #############
lambda0 <- 0.001 # penalization coefficient
psistars <- antitonicScoreMatch(Residualsorted, lambda0)
psistarprimes <- psistars[["psistarprimes"]]
psistars <- psistars[["psistars"]]

deltaknots <- diff(Residualsorted)
betahat <-spline_linear_regression_newton(betapilot, X, Y, Residualsorted, psistars, psistarprimes, 
                                          n_iter = 8, robust = robust, deltaknots = deltaknots)

# print(betahat[-1]-beta_0)
cat("Spline based SCALE's squared error ", sum((betahat[-1]-beta_0)^2), '\n')
# print(betapilot[-1]-beta_0)
cat("OLS's squared error ", sum((betapilot[-1]-beta_0)^2), '\n')
# print(betahat[1]-mu)
# print(betapilot[1]-mu)

############# Method B with CV#############
lambda0grid = 0.1^(0:4)
robust = FALSE

X1 <- X[1:(N / 2), ]
Y1 <- Y[1:(N / 2)]
X2 <- X[(N / 2 + 1):N, ]
Y2 <- Y[(N / 2 + 1):N]

X2temp <- cbind(rep(1, N/2), X2)
betapilot2 <- solve(t(X2temp) %*% X2temp, t(X2temp) %*% Y2) # betapilot2[1] is the intercept
Residualsorted1 <- sort(Y1 - X1 %*% betapilot2[-1]-betapilot2[1])

MSE <- rep(0, length(lambda0grid))
for (i in 1:length(lambda0grid)) {
  psistars <- antitonicScoreMatch(Residualsorted1, lambda0grid[i])
  psistarprimes <- psistars[["psistarprimes"]]
  psistars <- psistars[["psistars"]]
  
  deltaknots1 <- diff(Residualsorted1)
  betahat1 <-spline_linear_regression_newton(betapilot2, X1, Y1, Residualsorted1, psistars, psistarprimes, 
                                            n_iter = 8, robust = robust, deltaknots = deltaknots1)
  MSE[i] <- sum((X2 %*% (betahat1[-1]-betapilot2[-1]))^2)
}
cat("The lambda selected", lambda0grid[which.min(MSE)], "\n")
lambda0 <- lambda0grid[which.min(MSE)]
psistars <- antitonicScoreMatch(Residualsorted, lambda0)
psistarprimes <- psistars[["psistarprimes"]]
psistars <- psistars[["psistars"]]

deltaknots <- diff(Residualsorted)
betahat <-spline_linear_regression_newton(betapilot, X, Y, Residualsorted, psistars, psistarprimes, 
                                          n_iter = 8, robust = robust, deltaknots = deltaknots)

# print(betahat[-1]-beta_0)
cat("Spline based SCALE's error ", sum((betahat[-1]-beta_0)^2)^0.5, '\n')
# print(betapilot[-1]-beta_0)
cat("OLS's error ", sum((betapilot[-1]-beta_0)^2)^0.5, '\n')

