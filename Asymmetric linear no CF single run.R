source("scoreFun.R")

Noise_distr <- "uniform"
N <- 100# sample size
d <- 2 # dimension of X
beta_0 <- 3 * runif(d) + 1
mu <- 2

## Generate the data
X <- rnorm(N * d, 0, 1)
X <- array(X, c(N, d))

if(Noise_distr=="Gaussian"){Noise <- rnorm(N, 0, 1)}
if(Noise_distr=="uniform"){Noise <- runif(N, -1, 1)}
if(Noise_distr=="Cauchy"){Noise <- rcauchy(N)}
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

# ##use OLS as pilot
# betapilot <- solve(t(X_intercept) %*% X_intercept, t(X_intercept) %*% Y) # betapilot[1] is the intercept
# # betapilot[1] <- betapilot[1] + 3

###use median regressor as pilot
beta_init <- c(mu, beta_0)
result <- optim(par = beta_init, fn = mae_loss, X = X_intercept, Y = Y, method = "L-BFGS-B")
betapilot <- result$par

Residualsorted <- sort(Y - X_intercept %*% betapilot)

############# Method A #############
kernel <- "gaussian"
truncation_lim <- NULL
set_to_zero <- FALSE

fcn_psi <- kde_decr_score_est(Residualsorted, k=1000, kernel=kernel,
                              kernel_pts=2^20, truncation_lim=truncation_lim,
                              set_to_zero=set_to_zero)

fcn_psi_deriv <- fcn_psi[["psi_deriv"]]
fcn_psi <- fcn_psi[["psi"]]

# plot(fcn_psi, xlim = c(-0.5,0.5))
# plot(fcn_psi_deriv, xlim = c(-0.6,0.8))

# X_intercept <- cbind(rep(1, N), X)
betahat <- betapilot
# betahat <- solve(t(X_intercept) %*% X_intercept, t(X_intercept) %*% Y)

residuals <- Y - X_intercept %*% betahat
psibeta <- t(X_intercept) %*% fcn_psi(residuals)
psiprimebeta <- fcn_psi_deriv(residuals)
psiprimebeta <- -t(X_intercept) %*% diag(psiprimebeta) %*% X_intercept

n_iter = 40
for (l in 1:n_iter) {
  alpha <- 1
  if(sum(psibeta^2)^0.5 < (d / N)^0.5 * 1e-9) {
    break
  }
  print("norm of the z equation")
  print(sum(psibeta^2)^0.5)
  # residuals <- Y - X_intercept %*% betahat
  # psibeta <- t(X_intercept) %*% fcn_psi(residuals)
  # print(sum(psibeta^2)^0.5)
  # psiprimebeta <- fcn_psi_deriv(residuals)
  # psiprimebeta <- -t(X_intercept) %*% diag(psiprimebeta) %*% X_intercept
  # print("condition number")
  # print(kappa(psiprimebeta))
  # if(kappa(psiprimebeta) < 20){
  #   Update <- solve(psiprimebeta, psibeta)
  # }else{
  #   Update <- solve(psiprimebeta + diag(d+1), psibeta)
  # }
  Update <- solve(psiprimebeta + diag(d+1), psibeta)
  print("norm of the update")
  print(sum(Update^2)^0.5)
  # line search
  betatemp <- betahat - alpha * Update
  residuals <- Y - X_intercept %*% betatemp
  psibetatemp <- t(X_intercept) %*% fcn_psi(residuals)
  while (sum(psibetatemp^2) > 1 * sum(psibeta^2)){
    print("shrink the step size")
    alpha <- 0.8 * alpha
    betatemp <- betahat - alpha * Update
    residuals <- Y - X_intercept %*% betatemp
    psibetatemp <- t(X_intercept) %*% fcn_psi(residuals)
  }

  betahat <- betatemp
  print("betahat")
  print(betahat)
  print("new z norm")
  print(sum(psibetatemp^2)^0.5)
  print("**************************************")
  psiprimebeta <- fcn_psi_deriv(residuals)
  psiprimebeta <- -t(X_intercept) %*% diag(psiprimebeta) %*% X_intercept
  psibeta <- psibetatemp
}
###check if betahat solve the Z-equation
residuals <- Y - X_intercept %*% betahat
psibeta <- t(X_intercept) %*% fcn_psi(residuals)
cat("After",l,"iterations,","psibetahat has Euclidean norm ", sum(psibeta^2)^0.5, '\n')
cat("KDE based SCALE's squared error ", sum((betahat[-1]-beta_0)^2), '\n')
cat("Pilot's squared error ", sum((betapilot[-1]-beta_0)^2), '\n')

##### for inference ##### 
info <- mean((fcn_psi(residuals))^2)



############# Method B #############
robust = FALSE

lambda0 <- 1e-9 # penalization coefficient
psistars <- antitonicScoreMatch(Residualsorted, lambda0)
psistarprimes <- psistars[["psistarprimes"]]
psistars <- psistars[["psistars"]]

####reduce the memory use
non_redundent_index <- diff(psistars) != 0
temp <- c(non_redundent_index, FALSE)
non_redundent_index <- c(TRUE, non_redundent_index)
non_redundent_index <- non_redundent_index | temp
# print(N-sum(non_redundent_index))

# knots <- Residualsorted[non_redundent_index]
# psistars1<-psistars[non_redundent_index]
# psistarprimes1<-psistarprimes[non_redundent_index]
# range1 <- range(Residualsorted)
# xs <- seq(from = range1[1]-0.1, to = range1[2]+0.1, length.out = 1e5)
# y1<- eval_quad_Score(xs, Residualsorted, psistars, psistarprimes, 
#                      robust = FALSE, derivative = FALSE)
# y2<- eval_quad_Score(xs, knots, psistars1, psistarprimes1, 
#                      robust = FALSE, derivative = FALSE)
# print(sum(abs(y1-y2)))
# 
# y1<- eval_quad_Score(xs, Residualsorted, psistars, psistarprimes, 
#                      robust = FALSE, derivative = TRUE)
# y2<- eval_quad_Score(xs, knots, psistars1, psistarprimes1, 
#                      robust = FALSE, derivative = TRUE)
# print(sum(abs(y1[["psiprimexs"]]-y2[["psiprimexs"]])))


betapilot11<- solve(t(X_intercept) %*% X_intercept, t(X_intercept) %*% Y)
betahat <-spline_linear_regression_newton(betapilot11, X, Y, Residualsorted[non_redundent_index],
                                          psistars[non_redundent_index], psistarprimes[non_redundent_index],
                                          n_iter = 40, robust = robust)

# print(betahat[-1]-beta_0)
cat("Spline based SCALE's squared error ", sum((betahat[-1]-beta_0)^2), '\n')
# print(betapilot[-1]-beta_0)
cat("OLS's squared error ", sum((betapilot[-1]-beta_0)^2), '\n')



############# Method B with CV for lambda tuning#############
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
# 
# # print(betahat[-1]-beta_0)
# cat("Spline SCALE SE", sum((betahat[-1]-beta_0)^2), '\n')
# # print(betapilot[-1]-beta_0)
# cat("Pilot SE", sum((betapilot[-1]-beta_0)^2), '\n')

##### for inference ##### 
#info <- mean((fcn_psi(residuals))^2)

############# Oracle convex loss #############
###need to specify the noise pdf evaluated on a grid first
if(Noise_distr=="Gaussian"){beta_ora <- solve(t(X_intercept) %*% X_intercept, t(X_intercept) %*% Y)}
if(Noise_distr=="uniform"){
  result <- optim(par = betapilot, fn = max_loss, X = X_intercept, Y = Y, method = "BFGS")
  beta_ora <- result$par
  #check if the coverage rate is 1
  # range <- 2
  # print(sum(abs(X_intercept%*% beta_ora - Y) > (range / 2)))
}
if(Noise_distr=="Cauchy"){
  m <- min(1/(10*N), 1/2000)
  grids <- seq(m, 1-m, m)
  quantilegrid <- qcauchy(grids)
  J <- dcauchy(quantilegrid)
  init_scores <- diff(J) / m #psihat\circle Qhat
  res = isoProjInterpolate(init_scores, (quantilegrid[-1] + quantilegrid[-length(quantilegrid)]) / 2)
  ######## to be finished: z equation solver
  res[["iso_fn"]]
  res[["iso_fn_deriv"]]
  
}


############# One step #############
n <- floor(N/2)
X1_intercept <- X_intercept[1:n, ]
Y1 <- Y[1:n]
X2_intercept <- X_intercept[(n + 1):N, ]
Y2 <- Y[(n + 1):N]
residualsorted1 <- sort(Y1 - X1_intercept %*% betapilot)
residualsorted2 <- sort(Y2 - X2_intercept %*% betapilot)
kde1 <- density(residualsorted1, kernel = kernel, n = 2^21)
kde2 <- density(residualsorted2, kernel = kernel, n = 2^21)

### f'/f
temp1 <- diff(kde1$y)/diff(kde1$x)
temp11 <- (kde1$y[-1] + kde1$y[-length(kde1$y)])/2
index1 <- temp11 != 0
psi1 <- numeric(length = (length(temp1)))
psi1[index1] <- temp1[index1]/temp11[index1]
psi1 <- approxfun((kde1$x[-1] + kde1$x[-length(kde1$x)])/2, psi1, method = "linear",rule=2)
temp2 <- diff(kde2$y)/diff(kde2$x)
temp22 <- (kde2$y[-1] + kde2$y[-length(kde2$y)])/2
index2 <- temp22 != 0
psi2 <- numeric(length = (length(temp2)))
psi2[index2] <- temp2[index2]/temp22[index2]
psi2 <- approxfun((kde2$x[-1] + kde2$x[-length(kde2$x)])/2, psi2, method = "linear",rule=2)

### test if this more stable than just using f'/f
psi1 <- score_given_density(kde1$x, kde1$y, k=max(1000, 2*N))
psi2 <- score_given_density(kde2$x, kde2$y, k=max(1000, 2*N))

score1 <- -t(X1_intercept) %*% psi2(residualsorted1)
score2 <- -t(X2_intercept) %*% psi1(residualsorted2)

beta_one_step <- betapilot + solve(score1 %*% t(score1) + score2 %*% t(score2), score1 + score2)





