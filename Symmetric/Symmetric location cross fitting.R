library(ggplot2)
library(gridExtra)
library(gnorm) #to generate generalized Gaussian

# rm(list = ls())  # clear all variables in workspace
# graphics.off()  # clear all plots

N <- 1000 # sample size of each fold
K <- 3 # number of folds
lambda0 <- 1e-4 # penalization coefficient

set.seed(1)
## Generate data
X <- rexp(N*K)
X <- X * sign(runif(N, -1, 1)) # Laplace

# X <- rnorm(K*N, 0, 1)
# X <- runif(N*K, -1, 1) # U[-1, 1]
# X <- rgnorm(N*K, mu = 0, alpha = 1, beta = 9) # generalized Gaussian
# X <- runif(N*K, 0.5, 1) + runif(N, 0.5, 1) # triangular density: log - concave; infinite Fisher information
# X <- runif(N*K, -1, 1) + rnorm(N*K, 0, 1) # U[-1, 1] * N(0, 1)

# normalized penalization coefficient
lambda <- lambda0*2*N*sd(X)^3 #sd(X)^3: to make lambda scale invariant

#### K-fold cross-fitting: using one fold to estimate psi and 4 to estimate location, and switching roles
X <- matrix(X, nrow = K)
####### G[i, ] = (g*'(0), g*'(T1),..., g*'(Tn)) obtained on the i-th fold #######
G <- matrix(0, nrow = K, ncol = N+1)
Xsorteds <- matrix(0, nrow = K, ncol = N)
deltaXs <- matrix(0, nrow = K, ncol = N)
for (i in 1:K) {
  thetapilot <- mean(X[-i, ]) # pilot estimate: mean of the observations not in the i-th fold
  Xsorteds[i, ] <- sort(abs(X[i, ] - thetapilot)) # sorted out-of-sample residuals
  Xsorted <- Xsorteds[i, ]
  deltaXs[i, ] <- diff(c(0, Xsorted))
  deltaX <- deltaXs[i, ]
  u <- c(0, rep(1, N))
  LD1 <- t(matrix(replicate(N, deltaX), ncol = N))
  LD1[upper.tri(LD1)]  <-  0
  A <- diag(x = 1, nrow = N, ncol = N+1)
  diag(A[, -1]) <- 1
  D2 <- diag(1/deltaX*lambda) # lambda = lambda0*2*N*sd(X)^3
  B <- diag(x = 1, nrow = N, ncol = N+1)
  diag(B[, -1]) <- -1
  Hess <- 0.25*crossprod(LD1%*%A) + t(B)%*%D2%*%B
  ## initialization
  Gtemp <- solve(Hess, u)
  Gtemp[Gtemp < 0] <- 0 # projection
  M <- 15 # max number of iterations
  ## Projected Newton (main loop)
  for (j in 1:M) {
    # cat("Newton loop loss", t(Gtemp)%*%Hess%*%Gtemp - 2*t(u)%*%Gtemp,"i = ", i, "\n")
    alpha <- 1 # reset step size
    Free <- (Gtemp > 0) | (Hess%*%Gtemp - u <= 0)
    ### backtracking
    Gtemp1 <- Gtemp
    gradient <- (Hess%*%Gtemp - u)[Free]
    Update <- solve(Hess[Free,Free], gradient)
    Gtemp1[Free] <- Gtemp[Free] - alpha*Update
    Gtemp1[Free][Gtemp1[Free] < 0] <- 0 # projection
    while (t(Gtemp1)%*%Hess%*%Gtemp1 - 2*t(u)%*%Gtemp1 > t(Gtemp)%*%Hess%*%Gtemp - 2*t(u)%*%Gtemp + 0.5*t((Gtemp1-Gtemp)[Free])%*%(2*gradient)) {
      alpha <- 0.7*alpha
      Gtemp1[Free] <- Gtemp[Free] - alpha*Update
      Gtemp1[Free][Gtemp1[Free] < 0] <- 0 # projection
    }
    #### stopping criterion
    if (sum((Gtemp-Gtemp1)^2) < N*1e-12){
      Gtemp <- Gtemp1
      break
    }
    Gtemp <- Gtemp1
  }
  G[i, ] <- Gtemp
}

####### gvalues[i, ] = (gi*(0) = 0, gi*(T1),..., gi*(Tn)) obtained on the i-th fold #######
# K score function estimates, one (row) for each fold
gvalues <- G[, -1]
gvalues <- gvalues + G[, -(N+1)]
gvalues <- 0.5*gvalues*deltaXs
for (j in 2:N) {
  gvalues[, j] <- gvalues[, j] + gvalues[, j-1]
}
gvalues <- cbind(rep(0, K), gvalues)

######## dG_k/dx  = gk*, Negloglikelihd[k, ] = (Gk(0), Gk(T1),..., Gk(Tn))
Negloglikelihd <- gvalues[, -(N+1)]*deltaXs + (G[, -(N+1)]/2 + (G[, -1]-G[, -(N+1)])/6)*deltaXs^2
for (j in 2:N) Negloglikelihd[, j] <- Negloglikelihd[, j] + Negloglikelihd[, j-1]
normalization <- rep(0, K) # normalization[i] = log(integral exp(-G)) where G(0) = 0
Negloglikelihd <- sweep(Negloglikelihd, 1, normalization, "+")
Negloglikelihd <- cbind(normalization, Negloglikelihd)

############# PLOT ############
for (plotindex in 1:K) {
  tailgrid <- (1:(N/100))*mean(deltaXs[plotindex, ])
  Xplot <- c(0, Xsorteds[plotindex, ], tailgrid + Xsorteds[plotindex, N])
  ############# g ############
  NStailY <- gvalues[plotindex, N+1] + tailgrid*G[plotindex, N+1]
  NSYplot <- c(gvalues[plotindex, ], NStailY)
  ############# negative log-likelihood ############
  NLLtailY <- Negloglikelihd[plotindex, N+1] + tailgrid*gvalues[plotindex, N+1] + tailgrid^2/2*G[plotindex, N+1]
  NLLYplot <- c(Negloglikelihd[plotindex, ], NLLtailY)
  
  mydf <- data.frame(X = Xplot, NScore = NSYplot, NLL = NLLYplot)
  plotNS.obj <- ggplot(data = mydf) + geom_line(aes(x = X, y = NScore), color = "red") + ylab("Estimated negative score function") 
  # + geom_line(aes(x = X, y = 10*X^9), color = "blue") 
  plotNLL.obj <- ggplot(data = mydf) + geom_line(aes(x = X, y = NLL), color = "red") + ylab("Estimated negative log-likelihood")
  # + geom_line(aes(x = X, y = X^10), color = "blue")
  grid.arrange(plotNS.obj, plotNLL.obj, ncol = 2)
}
####################################

gtheta <- function(theta, X1, cuts, gvalues, G, deltaX){
  # deltaX <- diff(cuts)
  absdelta <- abs(X1 - theta)
  index <- findInterval(absdelta, cuts)
  Values <- gvalues[index] # g(left end point)
  distoleft <- absdelta - cuts[index]
  Values <- Values + distoleft*G[index] # linear part
  quadpts <- index != length(cuts)
  quadindex <- index[quadpts]
  Values[quadpts] <- Values[quadpts] + 0.5*(G[quadindex+1]-G[quadindex])*distoleft[quadpts]^2/deltaX[quadindex]
  Values <- Values*sign(X1 - theta)
  return(sum(Values))
}
### gprimetheta is dg/dx =-dg/dtheta
gprimetheta <- function(theta, X1, cuts, gvalues, G, deltaX){
  # deltaX <- diff(cuts)
  absdelta <- abs(X1 - theta)
  index <- findInterval(absdelta, cuts)
  Values <- G[index] # g(left end point)
  quadpts <- index != N+1
  quadindex <- index[quadpts]
  distoleft <- absdelta - cuts[index]
  Values[quadpts] <- Values[quadpts] + (G[quadindex+1]-G[quadindex])*distoleft[quadpts]/deltaX[quadindex]
  return(sum(Values))
}

#################### Cross-fitting: compute estimates thetahat_k for each fold and average ####################
thetahat1 <- 0
for (k in 1:K) {
  thetatemp <- 0 # initialization
  X1 <- c(X[-k, ])
  cuts <- c(0, Xsorteds[k, ])
  g <- gvalues[k, ]
  Gk <- G[k, ]
  deltaX <- deltaXs[k, ]
  ## Newton's method:
  while (abs(gtheta(thetatemp, X1, cuts, g, Gk, deltaX)) > 1e-13*N) {
    # #### backtracking ####
    # alpha <- 1
    # gvalue1 <- gtheta(thetatemp, X1, cuts, g, Gk, deltaX)
    # Update <- gvalue1/gprimetheta(thetatemp, X1, cuts, g, Gk, deltaX)
    # # gprimetheta is dg/dx = -dg/dtheta
    # thetatemp1 <- thetatemp + alpha*Update
    # while (abs(gtheta(thetatemp1, X1, cuts, g, Gk, deltaX)) > abs(gvalue1)) {
    #   alpha <- 0.8*alpha
    #   thetatemp1 <- thetatemp + alpha*Update
    # }
    # thetatemp <- thetatemp1
    #### no backtracking ####
    gvalue1 <- gtheta(thetatemp, X1, cuts, g, Gk, deltaX)
    Update <- gvalue1/gprimetheta(thetatemp, X1, cuts, g, Gk, deltaX)
    # gprimetheta is dg/dx =-dg/dtheta
    thetatemp <- thetatemp + Update
    if(abs(Update) < 1e-12) break
  }
  # print(abs(gtheta(thetatemp, X1, cuts, g, Gk, deltaX)))
  thetahat1 <- thetahat1 + thetatemp
}

thetahat1 <- thetahat1/K # averaged cross-fitted estimator
thetatilde <- mean(X) # compare with the sample mean

### Alternative to cross-fitting
# Use Newton's method to minimize the sum of the (estimated) loss function contributions from all K folds
#################### zero of sum_k g_k(X^{(-k)} - theta) ####################
#### no backtracking ####
thetahat2 <- 0
for (J in 1:12) {
  gthetahat <- 0
  gprimethetahat <- 0
  for (k in 1:K) {
    X1 <- c(X[-k, ])
    cuts <- c(0, Xsorteds[k, ])
    g <- gvalues[k, ]
    Gk <- G[k, ]
    deltaX <- deltaXs[k, ]
    gthetahat <- gthetahat + gtheta(thetahat2, X1, cuts, g, Gk, deltaX)
    gprimethetahat <- gprimethetahat + gprimetheta(thetahat2, X1, cuts, g, Gk, deltaX)
  }
  thetahat2 <- thetahat2 + gthetahat/gprimethetahat
  # cat("gthetahat = ", gthetahat,"\n")
  # cat("thetahat = ", thetahat,"\n")
}

cat("Sample mean =", thetatilde, ", thetahat_CF =", thetahat1, ", thetahat_Newton =", thetahat2)
################################################################################