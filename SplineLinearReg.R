## Score-aided estimation for linear regression based on two score estimation methods:
# 1. Monotone quadratic splines minimizing a penalized score matching objective
# 2. Antitonic projected cubic spline estimates
## No cross-fitting employed here
library(ggplot2)
library(gridExtra)
library(gnorm) # to generate generalized Gaussians
library(Matrix)
library(Iso)
library(stats)

# rm(list = ls())  # clear all variables in workspace
# graphics.off()  # clear all plots

N <- 7 # sample size of each fold
d <- 2 # dimension (p) of the design matrix X
beta_0 <- 3 * runif(d) + 1 # regression coefficients
mu <- 2 # intercept (population mean of Y)

lambda_0 <- 0.0001 # penalization coefficient

### Generate data
## Design matrix:
X <- matrix(rnorm(N * d, 0, 1), nrow = N)
# X <- array(X, c(N, d))
## normalized penalization coefficient
lambda <- lambda_0 * 2 * N * sd(X)^3 # sd(X)^3: to make lambda scale invariant

## Additive noise:
# Noise <- rcauchy(N)
# Noise <- runif(N, -1, 1)
Noise <- rnorm(N, 0, 1)
## log-concave, mean = median = 0, asymmetric
# Noise <- rnorm(N, 0, 1)
# negNoise <- (Noise < 0)
# Noise <- runif(N, - 2 * sqrt(2/pi), 0) * negNoise + Noise * (1-negNoise)
Y <- X %*% beta_0 + mu + Noise

### Step 1: Pilot estimation using OLS
Xtemp <- cbind(rep(1, N), X)
beta_pilot <- solve(t(Xtemp) %*% Xtemp, t(Xtemp) %*% Y) #beta_pilot[1] is the intercept
Residual_sorted <- sort(Y - Xtemp %*% beta_pilot)
delta_Residual <- diff(Residual_sorted)

### Step 2: Now use the residuals from Step 1 to estimate the score function
############ Shape-constrained splines
## Initialization
u <- c(0, rep(1, N))
L <- t(matrix(replicate(N, c(2, delta_Residual)), ncol = N))
L[upper.tri(L)] <- 0
A <- diag(x = 1, nrow = N, ncol = N+1)
diag(A[-1, 3:(N + 1)]) <- 1
A[1, 1] <- -1
D <- diag(1 / delta_Residual * lambda) #lambda =  lambda_0 * 2 * N * sd(X)^3
B <- diag(x = 1, nrow = N-1, ncol = N)
diag(B[, -1]) <- -1
#D and B above are not the B and D in the tex file;
#we only use them to calculate the non-zero lower right block of the Hessian matrix
Hess <- t(B) %*% D %*% B #non-zero lower right block of the Hessian matrix
Hess <- cbind(rep(0, N),Hess)
Hess <- rbind(rep(0, N + 1),Hess)
Hess <- 0.25 * crossprod(L %*% A) + Hess 
# Here, Hess = 0.5 * Hessian of the objective function, (Hess %*% G - u) = 0.5 * gradient
# G = (-g(T_1), g'(T_1), g'(T_2),..., g'(T_N))
# g = -psi is monotone increasing
G <- solve(Hess, u)
G[G < 0] <- 0 # projection

## Projected Newton (main loop):
M <- 9 # max number of iterations
for (j in 1:M) {
  # cat("Newton loop loss", t(G) %*% Hess %*% G - 2 * t(u) %*% G, "i = ", i,"\n")
  alpha <- 1 # reset step size
  Free <- (G > 0) | (Hess %*% G - u <= 0)
  ### backtracking
  Gtemp <- G
  gradient <- (Hess %*% G - u)[Free]
  Update <- solve(Hess[Free, Free], gradient)
  Gtemp[Free] <- G[Free] - alpha*Update
  Gtemp[Free][Gtemp[Free] < 0] <- 0 # projection
  while (t(Gtemp) %*% Hess %*% Gtemp - 2 * t(u) %*% Gtemp > t(G) %*% Hess %*% G - 2 * t(u) %*% G + 0.5 * t((Gtemp-G)[Free]) %*% (2*gradient)) {
    alpha <- 0.7 * alpha
    Gtemp[Free] <- G[Free] - alpha * Update
    Gtemp[Free][Gtemp[Free] < 0] <- 0 # projection
  }
  #### Termination criterion
  if (sum((G-Gtemp)^2) < N * 1e-13) {
    G <- Gtemp
    break
  }
  G <- Gtemp
}
####### G = (-g*(T1), g*'(T1),..., g*'(Tn)) #######

####### gvalues = (g*(T1),..., g*(Tn))#######
gvalues <- G[-1]
gvalues <- gvalues[-1] + gvalues[-N]
gvalues <- 0.5*gvalues*delta_Residual
for (j in 2:(N-1)) gvalues[j] <- gvalues[j] + gvalues[j-1]
gvalues <- c(-G[1], -G[1] + gvalues)
# sum(gvalues) # check == 0

####### gprimes = (g*'(T1),..., g*'(Tn)) #######
gprimes <- G[-1]

# Evaluate a quadratic monotone spline at a point x based on gvalues, gprimes, delta_Residual
gvalue <- function(x, knots, gvalues, gprimes, delta_Residual) {
  # delta_Residual <- diff(knots)
  index <- findInterval(x, knots)
  index[index == 0] <- 1 # left linear region
  Values <- gvalues[index] # g(left end point)
  distoleft <- x - knots[index]
  Values <- Values + distoleft*gprimes[index] # linear part
  quadpts <- (index != length(knots)) & (index != 1)
  quadindex <- index[quadpts]
  Values[quadpts] <- Values[quadpts] + 0.5 * (gprimes[quadindex + 1] - gprimes[quadindex]) * distoleft[quadpts]^2 / delta_Residual[quadindex]
  return(Values)
}

### Step 3: Computing the minimizer betahat of the convex loss function based on the estimated scores
## First derivative of the empirical risk: gbeta = - \nabla_{beta} sum_i loglikelihd(Y_i - <X_i,beta>), a vector
gbeta <- function(betahat, X, Y, knots, gvalues, gprimes, delta_Residual) {
  # delta_Residual <- diff(knots)
  X <- cbind(rep(1, length(Y)), X)
  residuals <- Y -X %*% betahat
  index <- findInterval(residuals, knots)
  index[index == 0] <- 1 # left linear region
  Values <- gvalues[index] # g(left end point)
  distoleft <- residuals - knots[index]
  Values <- Values + distoleft*gprimes[index] #linear part
  quadpts <- (index != length(knots)) & (index != 1)
  quadindex <- index[quadpts]
  Values[quadpts] <- Values[quadpts] + 0.5 * (gprimes[quadindex + 1] - gprimes[quadindex]) * distoleft[quadpts]^2 / delta_Residual[quadindex]
  # Values <- t(Values) %*% (X)
  # return(t(Values))
  return(-t(X) %*% Values)
}
## Second derivative of the empirical risk: gprimebeta = Hessian{ -sum_i loglikelihd(Y_i-<X_i,beta>)}
gprimebeta <- function(betahat, X, Y, knots, gvalues, gprimes, delta_Residual) {
  # delta_Residual <- diff(knots)
  X <- cbind(rep(1, length(Y)), X)
  delta <- (Y - X %*% betahat)
  index <- findInterval(delta, knots)
  index[index == 0] <- 1 # left linear region
  Values <- gprimes[index] # g(left end point)
  quadpts <- (index != length(knots)) & (index != 1)
  quadindex <- index[quadpts]
  distoleft <- delta - knots[index]
  Values[quadpts] <- Values[quadpts] + (gprimes[quadindex + 1] - gprimes[quadindex]) * distoleft[quadpts] / delta_Residual[quadindex]
  Values <- t(X) %*% diag(Values) %*% X
  return(Values)
}

## Compute the minimizer of the empirical risk using Newton's method
betahat <- beta_pilot # initialization: pilot OLS estimator from Step 1
for (l in 1:10) {
  gvalue1 <- gbeta(betahat, X, Y, Residual_sorted, gvalues, gprimes, delta_Residual)
  Update <- solve(gprimebeta(betahat, X, Y, Residual_sorted, gvalues, gprimes, delta_Residual), gvalue1)
  if (sum(abs(Update)) < d * 1e-9){
    betahat <- betahat - Update
    # print(sprintf("converged at %.0f th iteration", l))
    break
  }
  betahat <- betahat - Update
}
# print(sum(gbeta(betahat, X, Y, Residual_sorted, gvalues, G[-1], delta_Residual)^2))


########## Transforming an initial score estimate #############
# See Yu-Chun's spline note for details of the numerical implementation below
## First compute an initial cubic spline estimate (not monotone in general) 
A <- diag(N)
diag(A[-1, ]) <- -1

B <- diag(N)
diag(B[-1, ]) <- 1

W <- diag(c(0,1/delta_Residual))

v <- N:1
AAt_inv <- pmin(matrix(v, nrow = N, ncol = length(v), byrow = TRUE), matrix(v, nrow = N, ncol = length(v), byrow = FALSE)) # solve(A%*%t(A))
A_Binv <- matrix(2, ncol = N, nrow  = N)
A_Binv[upper.tri(A_Binv, diag = TRUE)] <- 0
A_Binv <- A_Binv + diag(N)
A_Binv <- A_Binv*sign(toeplitz(1:N%%2) - 0.5) #A %*% solve(B)

TEMP <- diag(c(1, delta_Residual^2)) %*% AAt_inv + 12 * lambda * W #=:[[r,v_1^t],[v_2, M]] see spline_opt.tex
Schur <- TEMP[-1, -1]#=M
Schur <- solve(Schur-matrix(TEMP[-1, 1]) %*% t(matrix(TEMP[1, -1]))/TEMP[1, 1])#=S see spline_opt.tex
TEMP <- -t(matrix(TEMP[1, -1])) %*% Schur / TEMP[1,1]
TEMP <- 6 * lambda * rbind(TEMP, Schur)
TEMP <- cbind(rep(0, N), TEMP)#=H_{11}^{-1}H_{22} see spline_opt.tex
# gstars = (g*(T1),..., g*(Tn)) is the negative score
# gstarprimes = (g*'(T1),..., g*'(Tn))
gstarprimes <- solve(lambda * (t(A_Binv) %*% W %*% A_Binv + 3 * W) - 6 * lambda * W^2 %*% TEMP, (N:1)%%2)
gstars <- TEMP %*% gstarprimes
gstarprimes <- cumsum(gstarprimes * sign((1:N)%%2 - 0.5)) * sign((1:N)%%2 - 0.5) #solve(B) %*% gstarprimes
gstars <- cumsum(gstars) # solve(A) %*% gstars
# sum(gstars) # check == 0

# Evaluate a cubic spline at a point x based on vectors knots, gstars, gstarprimes
gstar <- function(x, knots, gstars, gstarprimes) {
  N <- length(knots)
  deltaknots <- diff(knots)
  secantslopes <- diff(gstars) / deltaknots
  sumgprime <- gstarprimes[-N] + gstarprimes[-1]
  sumgprime1 <- 2 * gstarprimes[-N] + gstarprimes[-1]
  gvalues <- rep(0, length(x))
  index <- findInterval(x, knots)
  gvalues[index == 0] <- gstars[1] + gstarprimes[1] * (x[index == 0] - knots[1])
  gvalues[index == N] <- gstars[N] + gstarprimes[N] * (x[index == N] - knots[N])
  cubicpts <- !(index %in% c(0, N))
  cubicindex <- index[cubicpts]
  y <- (x[cubicpts] - knots[cubicindex]) / deltaknots[cubicindex]
  temp <- gstarprimes[cubicindex] * y + (3 * secantslopes[cubicindex] - sumgprime1[cubicindex]) * y^2 + (sumgprime[cubicindex] - 2 * secantslopes[cubicindex]) * y^3
  temp <- temp * deltaknots[cubicindex] + gstars[cubicindex]
  gvalues[cubicpts] <- temp
  return(gvalues)
}

# Next 'monotonize' the initial score function via isotonic regression (PAVA)
altindex <- (-1)^(1:((2 * N) - 1)) < 0
grid <- rep(0, 2 * N - 1)
grid[altindex] <- Residual_sorted
grid[!altindex] <- Residual_sorted[-N] + delta_Residual / 2
gIsos <- gstar(grid, Residual_sorted, gstars, gstarprimes)
gIsos <- pava(gIsos, decreasing = FALSE, long.out = FALSE, stepfun = FALSE)
gIso <- approxfun(grid, gIsos, rule = 2) # defines a function that can be evaluated

# gIso <- function(x, knots, gIsos){
#   gvalues <- rep(0, length(x))
#   index <- findInterval(x, knots)
#   gvalues[index == 0] <- gIsos[1]
#   gvalues[index] <- gIsos[index] #r ignores index 0
#   gvalues
# }

# gIsobeta <- function(betahat, X, Y, knots, gIsos){
#   X <- cbind(rep(1, length(Y)), X)
#   residuals <- Y - X %*% betahat
#   gvalues <- rep(0, length(residuals))
#   index <- findInterval(residuals, knots)
#   gvalues[index == 0] <- gIsos[1]
#   gvalues[index] <- gIsos[index]#r ignores index 0
#   -t(X)%*%gvalues
# }

## First derivative of the empirical risk
gIsobeta <- function(betahat, X, Y){
  X <- cbind(rep(1, length(Y)), X)
  residuals <- Y - X %*% betahat
  gvalues <- gIso(residuals) # evaluate the isotonic score estimate gIso computed above
  -t(X) %*% gvalues
}

secantslopes <- diff(gIsos) / diff(grid)
gIsoprime <- approxfun(grid, c(secantslopes, 0), method = "constant", yleft = 0, yright = 0, f = 0, ties = mean)

## Second derivative of the empirical risk
gIsobetaprime <- function(betahat, X, Y){
  X <- cbind(rep(1, length(Y)), X)
  residuals <- Y - X %*% betahat
  gprimes <- gIsoprime(residuals) # evaluate the derivative of the function gIso computed above
  return(t(X) %*% diag(gprimes) %*% X)
}

## Compute the minimizer of the empirical risk using Newton's method
betaIsohat <- beta_pilot # Initialization
for (l in 1:9) {
  Update <- solve(gIsobetaprime(betaIsohat, X, Y), gIsobeta(betaIsohat, X, Y))
  if (sum(abs(Update)) < d * 1e-9){
    betahat <- betahat - Update
    #print(sprintf("converged at %.0f th iteration", l))
    break
  }
  betaIsohat <- betaIsohat - Update
}
print(norm(gIsobeta(betaIsohat, X, Y)))

## Compare the estimation error of the three different methods
t(abs(beta_pilot - c(mu, beta_0))) # OLS
t(abs(betahat - c(mu, beta_0))) # Quadratic monotone spline
t(abs(betaIsohat - c(mu, beta_0))) # Projected spline estimate

psi_cauchy <- function(x){
  if (length(x) > 1) return(sapply(x, psi_cauchy))
  t0 <- newton(function(t) t * sin(t) + cos(t) - 1, 2)$root
  x0 <- cot(t0 / 2)
  x <- c(x, max(min(x, x0), -x0))
  # Return psi_0(x) and psi_0^*(x)
  -2*x / (1 + x^2)
}

############# Plot estimated score functions ############
# xgrid <- (-1800:3200)/1000 # grid of x values
xgrid <- seq(-4, 4, length.out = 5000)
Xplot <- c(Residual_sorted, xgrid)
# Score <- (sign(X) + 1)*X/2
# Score <- -psi_cauchy(Xplot)
Score <- Xplot
SplineY <- c(gstars, gstar(xgrid, Residual_sorted, gstars, gstarprimes))
SplineIsoY <- c(gIsos[altindex], gIso(xgrid))
SCSplineY <- c(gvalues, gvalue(xgrid, Residual_sorted, gvalues, gprimes, delta_Residual))
mydf <- data.frame(X = Xplot, Score = Score, SplineY = SplineY, SplineIsoY = SplineIsoY, SCSplineY = SCSplineY)
plot0.obj <- ggplot(data = mydf) + xlab("x") + ylab("y") +
  geom_line(aes(x = X, y = -Score, color = "Score")) +
  geom_line(aes(x = X, y = -SplineY, color = "SplineY")) + 
  geom_line(aes(x = X, y = -SplineIsoY, color = "SplineIsoY")) + 
  geom_line(aes(x = X, y = -SCSplineY, color = "SCSplineY")) +
  scale_color_manual(name = NULL,
                     breaks = c("Score", "SplineY", "SplineIsoY", "SCSplineY"),
                     values = c("red", "gold", "green", "blue"),
                     labels = c("Gaussian score function", "Cubic spline", "Projected cubic spline", "Monotone spline estimate")) +
  # geom_point(aes(x = Residual_sorted[1], y = gstars[1]), colour = "black") + 
  # geom_point(aes(x = Residual_sorted[N], y = gstars[N]), colour = "black")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom")
grid.arrange(plot0.obj)
####################################
