## Score-aided estimation for linear regression with symmetric errors, based on
# penalized score matching via monotone quadratic splines
## Two types of cross-fitting corresponding to DML1 and DM2 in Chernozhukov et al.
# betahat1: Averaging K estimates betahat_k based on different fold assignments
# betahat2: Optimizing a single objective function obtained by averaging K different empirical risk functions / estimating equations
library(ggplot2)
library(gridExtra)
library(gnorm) # to generate generalized Gaussians

# rm(list = ls())  # clear all variables in workspace
# graphics.off()  # clear all plots

N <- 500 # sample size of each fold
K <- 3 # number of folds for cross-fitting
d <- 2 # dimension of X
beta_0 <- 3*runif(d) + 1 # regression coefficients
mu <- 2 # intercept

lambda0 <- 1e-5 # penalization coefficient

## Generate data
X <- rnorm(K * N * d, 0, 1)
X <- array(X, c(K, N, d)) # 3D array
Y <- matrix(0, nrow = N, ncol = K)
for (i in 1:K) Y[, i] <- X[, , i]%*%beta_0 + mu # for each of the folds
## Additive noise
Noise  <- runif(K*N, -1, 1) # U[-1, 1] noise
# Noise  <- rnorm(K*N, 0, 1) # N(0, 1) noise
Noise  <- matrix(Noise, nrow = N, ncol = K)
Y <- Y + Noise

#### K-fold cross-fitting: use one fold to estimate psi and 4 to estimate beta, and then switch roles
# Normalized penalization coefficient:
lambda <- lambda0*2*N*sd(X)^3 # sd(X)^3: to make lambda scale invariant

####### G[i, ] = (g*'(0), g*'(T1),..., g*'(TN)) obtained on the i-th fold #######
G <- matrix(0, nrow = K, ncol = N + 1)
Residualsorteds <- matrix(0, nrow = K, ncol = N) # sorted residuals
deltaResiduals <- matrix(0, nrow = K, ncol = N) # differences between sorted residuals
betapilots <- matrix(0, nrow = d+1, ncol = K)
for (i in 1:K) {
  Xtemp <- apply(X[, , -i], 2, identity) # row-bind each data sheet; no standardization of the data
  Xtemp <- cbind(rep(1, (K - 1) * N), Xtemp) # add an intercept column to get a (K-1)*N times (p+1) matrix
  Ytemp <- as.vector(Y[, -i])
  
  ### Step 1: Pilot estimation using OLS (based on data not (?) in the i-th fold)
  betapilots[, i] <- solve(t(Xtemp)%*%Xtemp, t(Xtemp)%*%Ytemp) # betapilot[1] is the intercept
  betapilot <- betapilots[, i] # betapilot[1] is the intercept estimate
  # Out-of-sample residuals (sorted)
  Residualsorteds[i, ] <- sort(abs(Y[, i] - X[, , i]%*%betapilot[-1] - betapilot[1]))
  Residualsorted <- Residualsorteds[i, ]
  deltaResiduals[i, ] <- diff(c(0, Residualsorted))
  deltaResidual <- deltaResiduals[i, ]
  
  ### Step 2: Now use the residuals from Step 1 to estimate the score function using monotone splines
  ## Initialization
  u <- c(0, rep(1, N))
  LD1 <- t(matrix(replicate(N, deltaResidual), ncol = N))
  LD1[upper.tri(LD1)] <- 0
  A <- diag(x = 1, nrow = N, ncol = N + 1)
  diag(A[, -1]) <- 1
  D2 <- diag(1/deltaResidual*lambda) # lambda <- lambda0*2*N*sd(X)^3
  B <- diag(x = 1, nrow = N, ncol = N + 1)
  diag(B[, -1]) <- -1
  Hess <- 0.25*crossprod(LD1 %*% A) + t(B) %*% D2 %*% B
  Gtemp <- solve(Hess, u)
  Gtemp[Gtemp < 0] <- 0 # projection
  
  M <- 15 # max number of iterations
  ## Projected Newton (main loop)
  for (j in 1:M) {
    # cat("Newton loop loss", t(Gtemp)%*%Hess%*%Gtemp - 2*t(u)%*%Gtemp,"i=", i,"\n")
    alpha <- 1 # reset step size
    Free <- (Gtemp > 0) | (Hess%*%Gtemp - u <= 0)
    ### backtracking
    Gtemp1 <- Gtemp
    gradient <- (Hess%*%Gtemp - u)[Free]
    Update <- solve(Hess[Free, Free], gradient)
    Gtemp1[Free] <- Gtemp[Free] - alpha*Update
    Gtemp1[Free][Gtemp1[Free] < 0] <- 0 # projection
    while (t(Gtemp1)%*%Hess%*%Gtemp1 - 2*t(u)%*%Gtemp1 > t(Gtemp)%*%Hess%*%Gtemp - 2*t(u)%*%Gtemp + 0.5*t((Gtemp1-Gtemp)[Free])%*%(2*gradient)) {
      alpha <- 0.7*alpha
      Gtemp1[Free] <- Gtemp[Free] - alpha*Update
      Gtemp1[Free][Gtemp1[Free] < 0] <- 0 # projection
    }
    ## Termination criterion
    if (sum((Gtemp-Gtemp1)^2) < N*1e-12){
      Gtemp <- Gtemp1
      break
    }
    Gtemp <- Gtemp1
  }
  G[i, ] <- Gtemp
}

####### gvalues[i, ] = (g*(0) = 0, g*(T1),..., g*(TN)) obtained on the i-th fold #######
# K score function estimates, one (row) for each fold
gvalues <- G[, -1]
gvalues <- gvalues + G[, -(N+1)]
gvalues <- 0.5*gvalues*deltaResiduals
for (j in 2:N) gvalues[, j] <- gvalues[, j] + gvalues[, j-1]
gvalues <- cbind(rep(0, K), gvalues)

############# Plot the estimated score functions ############
for (plotindex in 1:K) {
  tailgrid <- (1:(N/100))*mean(deltaResiduals[plotindex, ])
  tailY <- gvalues[plotindex, N + 1] + tailgrid*G[plotindex, N + 1]
  Xplot <- c(0, Residualsorteds[plotindex,], tailgrid + Residualsorteds[plotindex,N])
  Yplot <- c(gvalues[plotindex, ], tailY)
  mydf <- data.frame(X = Xplot, Y = Yplot)
  plot0.obj <- ggplot(data = mydf) + geom_point(aes(x = X, y = Y), color = "red") + geom_line(aes(x = X, y = 0), color = "blue")
  grid.arrange(plot0.obj)
}
####################################

### Step 3: Computing the minimizer betahat of the convex loss function based on the estimated scores
## First derivative of the empirical risk: gbeta = - \nabla_{beta} sum_i loglikelihd(Y_i-<X_i,beta>), a vector
gbeta <- function(betahat, X, Y, cuts, gvalues, G, deltaResidual){
  # deltaResidual <- diff(cuts)
  X <- cbind(rep(1,length(Y)), X)
  residual <- Y - X%*%betahat
  absdelta <- abs(residual)
  index <- findInterval(absdelta, cuts)
  Values <- gvalues[index] # g(left end point)
  distoleft <- absdelta - cuts[index]
  Values <- Values + distoleft*G[index] # linear part
  quadpts <- index != length(cuts)
  quadindex <- index[quadpts]
  Values[quadpts] <- Values[quadpts] + 0.5*(G[quadindex+1]-G[quadindex])*distoleft[quadpts]^2/deltaResidual[quadindex]
  Values <- Values*sign(residual)
  Values <- t(Values)%*%(X)
  return(t(Values))
}

## Second derivative of the empirical risk: gprimebeta = Hessian{ - sum_i loglikelihd(Y_i-<X_i,beta>)}
gprimebeta <- function(betahat, X, Y, cuts, gvalues, G, deltaResidual){
  #deltaResidual <- diff(cuts)
  X <- cbind(rep(1,length(Y)), X)
  absdelta <- abs(Y - X%*%betahat)
  index <- findInterval(absdelta, cuts)
  Values <- G[index] # g(left end point)
  quadpts <- index != N+1
  quadindex <- index[quadpts]
  distoleft <- absdelta - cuts[index]
  Values[quadpts] <- Values[quadpts] + (G[quadindex+1]-G[quadindex])*distoleft[quadpts]/deltaResidual[quadindex]
  Values <- -t(X)%*%diag(Values)%*%X
  return(Values)
}

#################### Cross-fitting: Computing estimates betahat_k for each fold and averaging ####################
Xtemp <- apply(X[, ,], 2, identity) # row-bind each data sheet; no standardization of the covariates (no subtraction of their sample mean)
Xtemp <- cbind(rep(1, K*N), Xtemp)
Ytemp <- as.vector(Y[ , ])
### Pilot OLS estimate based on all of the data
betaOLS <- solve(t(Xtemp)%*%Xtemp, t(Xtemp)%*%Ytemp)

betahat <- betaOLS # initialization
#### no backtracking ####
for (k in 1:K) {
  ####### gbeta <- function(betahat, X, Y, cuts, gvalues, G, deltaResidual)
  X1 <- apply(X[, , -k], 2, identity) # row bind each data sheet
  Y1 <- as.vector(Y[, -k])
  cuts <- c(0, Residualsorteds[k, ])
  g <- gvalues[k, ]
  Gk <- G[k, ]
  deltaResidual <- deltaResiduals[k, ]
  betatemp <- betaOLS # Pilot OLS estimator
  ## Compute the minimizer of the empirical risk using Newton's method
  for (l in 1:5) {
    alpha <- 1
    gvalue1 <- gbeta(betatemp, X1, Y1, cuts, g, Gk, deltaResidual)
    Update <- solve(gprimebeta(betatemp, X1, Y1, cuts, g, Gk, deltaResidual), gvalue1)
    betatemp <- betatemp - alpha*Update
  }
  # print(sum(gbeta(betatemp, X1, Y1, cuts, g, Gk, deltaResidual)^2))
  betahat <- betahat + betatemp
}
# Average the K different regression estimates (cross-fitting)
betahat <- betahat/K

### Alternative to cross-fitting
# Use Newton's method to minimize the sum of the (estimated) loss function contributions from all K folds
#################### zero of sum_k g_k(X^{(-k)} - theta) ####################
# Xtemp <- apply(X[, ,], 2, identity) # row-bind each data sheet; no standardization of the data
# Xtemp <- cbind(rep(1, K*N), Xtemp)
# Ytemp <- as.vector(Y[ , ])
# ### OLS as pilot estimator
# betaOLS <- solve(t(Xtemp)%*%Xtemp, t(Xtemp)%*%Ytemp)

#### no backtracking ####
betahat2 <- betaOLS # Initialization
for (J in 1:7) {
  gbetahat <- 0
  gprimebetahat <- rep(0, d+1)
  for (k in 1:K) {
    X1 <- apply(X[, , -k], 2, identity) # row-bind each data sheet
    Y1 <- as.vector(Y[, -k])
    cuts <- c(0, Residualsorteds[k, ])
    g <- gvalues[k, ]
    Gk <- G[k, ]
    deltaResidual <- deltaResiduals[k, ]
    gbetahat <- gbetahat + gbeta(betahat2, X1, Y1, cuts, g, Gk, deltaResidual)
    gprimebetahat <- gprimebetahat + gprimebeta(betahat2, X1, Y1, cuts, g, Gk, deltaResidual)
  }
  betahat2 <- betahat2 - solve(gprimebetahat, gbetahat)
  # cat("gbetahat = ", gbetahat,"\n")
  # cat("thetahat = ", thetahat,"\n")
  # print(sum(gbetahat^2))
}
abs(betahat - c(mu, beta_0))
abs(betahat2 - c(mu, beta_0))
abs(betaOLS - c(mu, beta_0))
################################################################################