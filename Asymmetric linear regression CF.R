## Score-aided estimation for linear regression based on
# penalized score matching via monotone quadratic splines
## Two types of cross-fitting corresponding to DML1 and DM2 in Chernozhukov et al.
# betahat1: Averaging K estimates betahat_k based on different fold assignments
# betahat2: Optimizing a single objective function obtained by averaging K different empirical risk functions / estimating equations
library(ggplot2)
library(gridExtra)
library(gnorm) #to generate generalized Gaussian

# rm(list = ls())  # clear all variables in workspace
# graphics.off()  # clear all plots

N <- 800 # sample size of each fold
K <- 3 # number of folds
d <- 4 # dimension of X
beta_0 <- 3*runif(d) + 1 
mu <- 2 # intercept
lambda0<-0.0001 # penalization coefficient

## Generate the data
X <- rnorm(K*N*d, 0, 1)
X <- array(X, c(N, d, K))
Y <- matrix(0, nrow = N, ncol = K)
for (i in 1:K) Y[, i] <- X[, , i]%*%beta_0 + mu
# Additive noise
Noise <- runif(K*N, -1, 1)
# Noise <- rnorm(K*N, 0, 1)
# Asymmetric log - concave noise with mean = median = 0
# Noise <- rnorm(K*N, 0, 1)
# negNoise <- Noise < 0
# Noise <- runif(K*N,  - 2 *sqrt(2/pi), 0)*negNoise + Noise*(1-negNoise)
Y <- Y + Noise

#### K-fold cross-fitting: using one fold to estimate psi and 4 to estimate beta, and switching roles
# normalized penalization coefficient
lambda <- lambda0*2*N*sd(X)^3 # sd(X)^3: to make lambda scale invariant
####### G[i, ] = (-g*(T1), g*'(T1),..., g*'(Tn)) obtained on the i-th fold #######
G <- matrix(0, nrow = K, ncol = N+1)
Residualsorteds <- matrix(0, nrow = K, ncol = N)
deltaResiduals <- matrix(0, nrow = K, ncol = (N-1))
betapilots <- matrix(0, nrow = d+1, ncol = K)
for (i in 1:K) {
  Xtemp <- apply(X[, , -i], 2, identity) # row-bind each data sheet; no standardization of the covariates (no subtraction of covariates' sample mean)
  Xtemp <- cbind(rep(1, (K-1)*N), Xtemp)
  Ytemp <- as.vector(Y[, -i])
  ### Pilot OLS estimator
  betapilots[, i] <- solve(t(Xtemp)%*%Xtemp, t(Xtemp)%*%Ytemp)# betapilot[1] is the intercept
  betapilot <- betapilots[, i]
  # Sorted out-of-sample residuals
  Residualsorteds[i, ] <- sort(Y[, i] - X[, , i]%*%betapilot[-1] - betapilot[1])
  Residualsorted <- Residualsorteds[i, ]
  deltaResiduals[i, ] <- diff(Residualsorted)
  deltaResidual <- deltaResiduals[i, ]
  u <- c(0, rep(1, N))
  L <- t(matrix(replicate(N, c(2, deltaResidual)), ncol = N))
  L[upper.tri(L)] <- 0
  A <- diag(x = 1, nrow = N, ncol = N+1)
  diag(A[-1, 3:(N+1)]) <- 1
  A[1, 1] <- -1
  D <- diag(1/deltaResidual*lambda) # lambda = lambda0*2*N*sd(X)^3
  B <- diag(x = 1, nrow = N-1, ncol = N)
  diag(B[, -1]) <- -1
  #D and B above are not the B and D in the tex file;
  #we only use them to calculate the non-zero lower right block of the Hessian matrix
  Hess <- t(B)%*%D%*%B #non-zero lower right block of the Hessian matrix
  Hess <- cbind(rep(0, N), Hess)
  Hess <- rbind(rep(0, N+1), Hess)
  Hess <- 0.25*crossprod(L%*%A) + Hess 
  # Hess = 0.5*Hessian of the objective function, (Hess%*%G - u) = 0.5*gradient
  ## Initialization
  # Gtemp = (-g(T_1), g'(T_1), g'(T_2),..., g'(T_N))
  # g = -psi is monotone increasing
  Gtemp <- solve(Hess, u)
  Gtemp[Gtemp < 0] <- 0 # projection
  ## Projected Newton (main loop):
  M <- 9 # max number of iterations
  for (j in 1:M) {
    # cat("Newton loop loss", t(Gtemp)%*%Hess%*%Gtemp - 2 *t(u)%*%Gtemp, "i=", i, "\n")
    alpha <- 1 # reset step size
    Free <- (Gtemp>0)|(Hess%*%Gtemp - u <= 0)
    ### backtracking
    Gtemp1 <- Gtemp
    gradient <- (Hess%*%Gtemp - u)[Free]
    Update <- solve(Hess[Free,Free], gradient)
    Gtemp1[Free] <- Gtemp[Free] - alpha*Update
    Gtemp1[Free][Gtemp1[Free] < 0] <- 0 # projection
    while (t(Gtemp1)%*%Hess%*%Gtemp1 - 2 *t(u)%*%Gtemp1 > t(Gtemp)%*%Hess%*%Gtemp - 2 *t(u)%*%Gtemp+0.5*t((Gtemp1-Gtemp)[Free])%*%(2*gradient)) {
      alpha <- 0.7*alpha
      Gtemp1[Free] <- Gtemp[Free] - alpha*Update
      Gtemp1[Free][Gtemp1[Free] < 0] <- 0 # projection
    }
    #### Termination criterion
    if (sum((Gtemp-Gtemp1)^2) < N*1e-13){
      Gtemp <- Gtemp1
      break
    }
    Gtemp <- Gtemp1
  }
  G[i, ] <- Gtemp
}

####### gvalues[i, ] = (g*(T1),..., g*(Tn)) obtained on the i-th fold#######
gvalues <- G[, -1]
gvalues <- gvalues[, -1] + gvalues[, -N]
gvalues <- 0.5*gvalues*deltaResiduals
for (j in 2:(N-1)) {
  gvalues[, j] <- gvalues[, j] + gvalues[, j-1]
}
gvalues <- sweep(gvalues, 1, -G[, 1], "+")
gvalues <- cbind(-G[, 1], gvalues)
# rowSums(gvalues)

############# Plot the estimated score functions g ############
for (plotindex in 1:K) {
  tailgrid <- (1:(N/100))*mean(deltaResiduals[plotindex, ])
  tailY <- gvalues[plotindex, N] + tailgrid*G[plotindex, N+1]
  Negtailgrid <- (-(N/100):-1)*mean(deltaResiduals[plotindex, ])
  NegtailY <- gvalues[plotindex, 1] + Negtailgrid*G[plotindex,2]
  Xplot <- c(Negtailgrid+Residualsorteds[plotindex, 1], Residualsorteds[plotindex, ], tailgrid+Residualsorteds[plotindex, N])
  Yplot <- c(NegtailY, gvalues[plotindex, ], tailY)
  
  mydf <- data.frame(X = Xplot, Y = Yplot)
  plot0.obj <- ggplot(data = mydf) + geom_point(aes(x = X, y = Y), color = "red") + geom_line(aes(x = X, y = (sign(X) + 1)*X/2), color = "blue")
  grid.arrange(plot0.obj)
}
####################################

### gbeta= - \nabla_{beta} sum_i loglikelihd(Y_i-<X_i,beta>), a vector
gbeta <- function(betahat, X, Y, cuts, gvalues, GG, deltaResidual){
  # deltaResidual <- diff(cuts)
  X <- cbind(rep(1, length(Y)), X)
  residual <- Y - X%*%betahat
  delta <- (residual)
  index <- findInterval(delta, cuts)
  index[index == 0] <- 1# left linear region
  Values <- gvalues[index] # g(left end point)
  distoleft <- delta - cuts[index]
  Values <- Values + distoleft*GG[index] # linear part
  quadpts <- (index != length(cuts)) & (index != 1)
  quadindex <- index[quadpts]
  Values[quadpts] <- Values[quadpts]+0.5*(GG[quadindex+1]-GG[quadindex])*distoleft[quadpts]^2/deltaResidual[quadindex]
  Values <- t(Values)%*%(X)
  return(t(Values))
}
### gprimebeta = Hessian{ - sum_i loglikelihd(Y_i-<X_i,beta>)}
gprimebeta <- function(betahat, X, Y, cuts, gvalues, GG, deltaResidual){
  # deltaResidual <- diff(cuts)
  X <- cbind(rep(1, length(Y)), X)
  delta <- (Y - X%*%betahat)
  index <- findInterval(delta, cuts)
  index[index == 0] <- 1# left linear region
  Values <- GG[index] # g(left end point)
  quadpts <- (index != length(cuts)) & (index != 1)
  quadindex <- index[quadpts]
  distoleft <- delta - cuts[index]
  Values[quadpts] <- Values[quadpts] + (GG[quadindex+1]-GG[quadindex])*distoleft[quadpts]/deltaResidual[quadindex]
  Values <- -t(X)%*%diag(Values)%*%X
  return(Values)
}
#################### Cross-fitting: compute estimates betahat_k for each fold and then average ####################
Xtemp <- apply(X[,, ], 2, identity) # row-bind each data sheet; no standardization of the data
Xtemp <- cbind(rep(1, K*N), Xtemp)
Ytemp <- as.vector(Y[, ])
### Pilot OLS estimator
betaOLS <- solve(t(Xtemp)%*%Xtemp, t(Xtemp)%*%Ytemp)

betahat <- rep(0, d) # Initialization
for (k in 1:K) {
  betatemp <- betaOLS # Initialization
  ####### gbeta <- function(betahat, X, Y, cuts, gvalues, G, deltaResidual)
  X1 <- apply(X[, , -k], 2, identity) # row-bind each data sheet
  Y1 <- as.vector(Y[, -k])
  cuts <- Residualsorteds[k, ]
  g <- gvalues[k, ]
  GGk <- G[k, -1]
  deltaResidual <- deltaResiduals[k, ]
  ## Newton's method:
  for (l in 1:10) {
    alpha <- 1
    gvalue1 <- gbeta(betatemp, X1, Y1, cuts, g, GGk, deltaResidual)
    Update <- solve(gprimebeta(betatemp, X1, Y1, cuts, g, GGk, deltaResidual), gvalue1)
    betatemp <- betatemp - alpha*Update
  }
  # print(sum(gbeta(betatemp, X1, Y1, cuts, g, GGk, deltaResidual)^2))
  betahat <- betahat + betatemp
}
betahat <- betahat/K

################################################################################
####################zero of sum_k g_k(X^{(-k)}-theta)####################
#### primitive, no backtracking####
# Xtemp <- apply(X[,, ], 2, identity)# row-bind each data sheet; no standardization of the data
# Xtemp <- cbind(rep(1, K*N), Xtemp)
# Ytemp <- as.vector(Y[, ])
# ###OLS as pilot estimator
# betaOLS <- solve(t(Xtemp)%*%Xtemp, t(Xtemp)%*%Ytemp)

betahat2 <- betaOLS # Initialization
for (J in 1:7) {
  gbetahat <- 0
  gprimebetahat <- rep(0, d+1)
  for (k in 1:K) {
    X1 <- apply(X[, , -k], 2, identity)# row-bind each data sheet
    Y1 <- as.vector(Y[, -k])
    cuts <- Residualsorteds[k, ]
    g <- gvalues[k, ]
    GGk <- G[k, -1]
    deltaResidual <- deltaResiduals[k, ]
    gbetahat <- gbetahat + gbeta(betahat2, X1, Y1, cuts, g, GGk, deltaResidual)
    gprimebetahat <- gprimebetahat + gprimebeta(betahat2, X1, Y1, cuts, g, GGk, deltaResidual)
  }
  betahat2 <- betahat2-solve(gprimebetahat, gbetahat)
  # cat("gbetahat2 = ", gbetahat2, "\n")
  print(sum(gbetahat^2))
}

sum(abs(betahat - c(mu,beta_0))^2)
sum(abs(betahat2 - c(mu,beta_0))^2)
sum(abs(betaOLS - c(mu,beta_0))^2)
################################################################################
