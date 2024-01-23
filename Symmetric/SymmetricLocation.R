#library(RSpectra) #svds() used in gradient descent
library(ggplot2)
library(gridExtra)
library(gnorm) #to generate generalized Gaussian
#library(quadprog) #solve.QP(), off-the-shelf quadratic programming

#rm(list = ls())  # clear all variables in workspace
#graphics.off()  # clear all plots

N <- 2000 # sample size
lambda0 <- 1e-5 # penalization coefficient

## Generate data:
# X <- rexp(N)
# X <- X*sign(runif(N,-1,1)) # Laplace

# X <- runif(N,-1,1) # U[-1,1]

X <- rnorm(N, 0, 1) # Gaussian

# X <- rgnorm(N, mu = 0, alpha = 1, beta = 9) # generalized Gaussian
# X <- runif(N, 0.5, 1) + runif(N, 0.5, 1) # triangular density: log-concave; infinite Fisher information
# X <- runif(N, -1, 1) + rnorm(N, 0, 1) # U[-1,1] * N(0,1)

## Add beta noise:
# Z <- rbeta(N, 5, 5) - 0.5
# X <- X + 0.05*Z
lambda <- lambda0*2*N*sd(X)^3 # sd(X)^3: to make lambda scale invariant

## Pilot estimate of the location parameter theta:
# thetapilot <- mean(X)
thetapilot <- 0 # oracle/true value
Xsorted <- sort(abs(X - thetapilot))
deltaX <- diff(c(0, Xsorted))

## small scale projected Newton
u <- c(0, rep(1, N)) # linear part of the objective function
LD1 <- t(matrix(replicate(N, deltaX), ncol = N))
LD1[upper.tri(LD1)] <- 0
A <- diag(x = 1, nrow = N, ncol = N+1)
diag(A[, -1]) <- 1
#### test if the matrix expression is correct ####
# G<-runif(N+1, min = 0.1, max = 0.2)
# tempp <- 0
# for (i in 1:N) {
#   temp2 <- 0
#   for (j in 1:i) {
#     temp2 <- temp2 + (G[j] + G[j+1])*deltaX[j]
#   }
#   tempp <- tempp + temp2^2
# }
# tempp - t(G)%*%crossprod(LD1%*%A)%*%G
####################################################

D2 <- diag(1/deltaX*lambda) #lambda = lambda0*2*N*sd(X)^3
B<-diag(x = 1, nrow = N, ncol = N+1)
diag(B[, -1]) <- -1

#### test if the matrix expression is correct ####
# G<-runif(N+1, min = 0.1, max = 0.2)
# tempp <- 0
# for (j in 1:N) {
#   tempp <- tempp+(G[j+1]-G[j])^2/deltaX[j]
# }
# tempp*lambda*N*2 - t(G)%*%t(B)%*%D2%*%B%*%G
####################################################

Hess <- 0.25*crossprod(LD1%*%A) + t(B)%*%D2%*%B 
# Here, Hess = 0.5*Hessian of the objective function, (Hess%*%G-u) = 0.5*gradient
# Objective function: t(G)%*%Hess%*%G - 2*t(u)%*%G
## Initialization, G = (g'(0), g'(T1),..., g'(Tn))
# g = -psi is monotone increasing
G <- solve(Hess, u)
G[G < 0] <- 0 # projection
### Projected Newton (main loop):
M <- 20
for (i in 1:M) {
  cat("Newton loop loss", t(G)%*%Hess%*%G - 2*t(u)%*%G, "\n")
  alpha <- 1 # reset step size
  # Fixed <- (G == 0) * (Hess%*%G - u > 0)
  # Free <- !Fixed
  Free <- (G > 0) | (Hess%*%G - u <= 0)

  ## no backtracking for alpha
  # G[Free]<-G[Free] - alpha*Scaling[Free,Free]%*%(Hess%*%G - u)[Free]
  # G[Free][G[Free] < 0] <- 0 #projection

  ## backtracking
  Gtemp <- G
  gradient <- (Hess%*%G - u)[Free]
  Update <- solve(Hess[Free, Free], gradient)
  Gtemp[Free] <- G[Free] - alpha*Update
  Gtemp[Free][Gtemp[Free] < 0] <- 0 # projection
  while (t(Gtemp)%*%Hess%*%Gtemp - 2*t(u)%*%Gtemp > t(G)%*%Hess%*%G - 2*t(u)%*%G + 0.5*t((Gtemp - G)[Free])%*%(2*gradient)) {
    alpha <- 0.8*alpha
    Gtemp[Free] <- G[Free] - alpha*Update
    Gtemp[Free][Gtemp[Free] < 0] <- 0 # projection
  }
  if (sum((Gtemp-G)^2) < 1e-30){
    cat("Last update has Euclidean norm", sum((Gtemp-G)^2))
    G <- Gtemp
    break
  }
  G <- Gtemp
}

# ## small scale projected gradient descent
# Stepsize <- 1/svds(Hess, 1, nu = 0, nv = 0)$d[1] # Hess = 0.5*Hessian
# L <- 40
# for (i in 1:L) {
#   if (i %% 10 == 0) {
#     cat("Gradient descent loop loss", t(G)%*%Hess%*%G - 2*t(u)%*%G, "\n")
#   }
#   G <- G - Stepsize*(Hess%*%G - u)
#   G[G < 0] <- 0 # projection
# }

########### off-the-shelf QP ######################
# same result
# qp <- solve.QP(Dmat = Hess, dvec = u, Amat = diag(x = 1, nrow = N+1, ncol = N+1), bvec = rep(0, N+1), meq = 0)
# sum((qp$solution-G)^2)
############################################

# Obtain the values of g from g'
gvalues <- G[-1] + G[-(N+1)]
gvalues <- 0.5*gvalues*deltaX
# Or use cumsum:
for (j in 2:N) {
  gvalues[j] <- gvalues[j] + gvalues[j-1]
}
gvalues <- c(0, gvalues) # gvalues = (g*(0), g*(T1),..., g*(Tn))

######## G' = g*, Negloglikelihd = (G(0), G(T1),..., G(Tn))
Negloglikelihd <- gvalues[-(N+1)]*deltaX + (G[-(N+1)]/2 + (G[-1]-G[-(N+1)])/6)*deltaX^2
for (j in 2:N) {
  Negloglikelihd[j] <- Negloglikelihd[j] + Negloglikelihd[j-1]
}
normalization <- 0 # normalization[i] = log(integral exp(-G)) where G(0) = 0
Negloglikelihd <- Negloglikelihd + normalization
Negloglikelihd <- c(0, Negloglikelihd)

#############PLOT############
tailgrid <- (1:(N/100))*mean(deltaX)
Xplot <- c(0, Xsorted, tailgrid + Xsorted[N])
############# g*: negative score function ############
NStailY <- gvalues[N+1] + tailgrid*G[N+1]
NSYplot <- c(gvalues, NStailY)
############# negative log-likelihood ############
NLLtailY <- Negloglikelihd[N+1] + tailgrid*gvalues[N+1] + tailgrid^2/2*G[N+1]
NLLYplot <- c(Negloglikelihd, NLLtailY)
Xmin <- min(Xsorted)
Xmax <- max(Xsorted)

mydf <- data.frame(X = Xplot, NScore = NSYplot, NLL = NLLYplot)
plotNS.obj <- ggplot(data = mydf) + geom_line(aes(x = X, y = NScore), color = "red") + geom_line(aes(x = X, y = X), color = "blue") + # true score function (blue) is only linear in the Gaussian case
  ylab("Estimated negative score function") + geom_vline(xintercept = Xmin, linetype = "dotted") + geom_vline(xintercept = Xmax, linetype = "dotted")

plotNLL.obj <- ggplot(data = mydf) + geom_line(aes(x = X, y = NLL), color = "red") + geom_line(aes(x = X, y = X^2/2), color = "blue") +
  ylab("Estimated negative log-likelihood") + geom_vline(xintercept = Xmin, linetype = "dotted") + geom_vline(xintercept = Xmax, linetype = "dotted")

grid.arrange(plotNS.obj, plotNLL.obj, ncol = 2)
####################################

# ####### gvalues = (g(0),g(T1),...,g(Tn)))) #######
# ####### G = (g'(0),g'(T1),...,g'(Tn))#######
# # cuts <- c(0,Xsorted)
# # gtheta <- function(theta,X1){
# #   absdelta <- abs(X1-theta)
# #   index <- findInterval(absdelta, cuts)
# #   Values <- gvalues[index] #g(left end point)
# #   distoleft <- absdelta-cuts[index]
# #   Values <- Values+distoleft*G[index] #linear part
# #   quadpts <- index!=N+1
# #   quadindex <- index[quadpts]
# #   Values[quadpts] <- Values[quadpts]+0.5*(G[quadindex+1]-G[quadindex])*distoleft[quadpts]^2/deltaX[quadindex]
# #   Values <- Values*sign(X1-theta)
# #   return(sum(Values))
# # }
# # 
#####gprimetheta is d/dx g=-d/dtheta g
# # gprimetheta <- function(theta,X1){
# #   absdelta <- abs(X1-theta)
# #   index <- findInterval(absdelta, cuts)
# #   Values <- G[index] #g(left end point)
# #   quadpts <- index!=N+1
# #   quadindex <- index[quadpts]
# #   distoleft <- absdelta-cuts[index]
# #   Values[quadpts] <- Values[quadpts]+(G[quadindex+1]-G[quadindex])*distoleft[quadpts]/deltaX[quadindex]
# #   return(sum(Values))
# # }
# 
# ####### Compare with the sample mean
# Thetahats <- NULL
# Thetatildes <- NULL
# for (j in 1:1000) {
#   #X1 <- runif(N,-1,1)
#   X1 <- rnorm(N)
#   #X1 <- rgnorm(N, mu = 0, alpha = 1, beta = 4)
#   thetatilde <- mean(X1)
#   thetahat <- thetatilde #initialization
# 
#   while (abs(gtheta(thetahat,X1))>1e-13*N) {
#     ###backtracking
#     # alpha <- 1
#     # Update <- gtheta(thetahat,X1)/gprimetheta(thetahat,X1)
#     # if(abs(Update)<1e-12){break}
#     # thetatemp <- thetahat+alpha*Update
#     # while (abs(gtheta(thetatemp,X1))>abs(gtheta(thetahat,X1))) {
#     #   alpha <- 0.8*alpha
#     #   thetatemp <- thetahat+alpha*Update
#     # }
#     # #print(alpha)
#     # #print(abs(gtheta(thetatemp,X1))-abs(gtheta(thetahat,X1)))
#     # thetahat <- thetatemp
#     ###no backtracking
#     Update <- gtheta(thetahat,X1)/gprimetheta(thetahat,X1)
#     thetahat <- thetahat+Update
#     if(abs(Update)<1e-12){break}
#   }
#   #print(abs(gtheta(thetahat,X1)))
#   Thetahats <- c(Thetahats,abs(thetahat))
#   Thetatildes <- c(Thetatildes,abs(thetatilde)/sqrt(2))
#   #divide the |mean-theta0| by sqrt(2) to balance the extra 1000 samples used for score matching
# }
# mean(Thetahats)
# mean(Thetatildes)
