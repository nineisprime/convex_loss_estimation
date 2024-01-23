## Spline-based score estimation plots based on iid observations from a symmetric density
# (not based on a location or regression model)
#library(RSpectra) # svds() used in gradient descent
library(ggplot2)
library(gridExtra)
library(pracma)
# library(gnorm) # to generate generalized Gaussian

# rm(list = ls())  # clear all variables in workspace
# graphics.off()  # clear all plots

# Estimation of a antitonic symmetric projected score function using the penalized spline method
## Main function for computing a quadratic monotone spline fit by Newton's method 
# (see Yu-Chun's monotone splines note)
## Return the derivative vector (g'(0),g'(T1),...,g'(Tn)), where g = -psi is increasing
symmetric_neg_decr_score_est <- function(X, theta_pilot = 0, lambda_0 = 1e-4) {
  N <- length(X)
  lambda <- lambda_0 * 2 * N * sd(X)^3 # sd(X)^3: to make lambda scale invariant
  X_sorted <- sort(abs(X - theta_pilot))
  delta_X <- diff(c(0, X_sorted))
  
  ## Small scale projected Newton:
  u <- c(0, rep(1, N)) # linear part of the objective function
  LD <- t(matrix(replicate(N, delta_X), ncol = N))
  LD[upper.tri(LD)] <- 0
  A <- diag(x = 1, nrow = N, ncol = N+1)
  diag(A[, -1]) <- 1
  
  Dinv <- diag(1/delta_X * lambda) # lambda = lambda_0 * 2 * N * sd(X)^3
  B <- diag(x = 1, nrow = N, ncol = N+1)
  diag(B[, -1]) <- -1
  
  Hess <- 0.25 * crossprod(LD %*% A) + t(B) %*% Dinv %*% B 
  # Here, Hess = 0.5 * Hessian of the objective function, (Hess %*% G - u) = 0.5 * gradient
  # Objective function: t(G) %*% Hess %*% G - 2 * t(u) %*% G
  ## Initialization, G = (g'(0), g'(T1),..., g'(Tn))
  # g = -psi is monotone increasing
  G <- solve(Hess, u)
  G[G < 0] <- 0 # projection onto the non-negative orthant
  
  ### Projected Newton main loop:
  M <- 20
  for (i in 1:M) {
    cat("Newton loop loss", t(G) %*% Hess %*% G - 2 * t(u) %*% G, "\n")
    alpha <- 1 # reset step size
    # Fixed <- (G == 0) * (Hess%*%G - u > 0)
    # Free <- !Fixed
    Free <- (G > 0) | (Hess %*% G - u <= 0)
    
    ### no backtracking for alpha
    # G[Free] <- G[Free] - alpha*Scaling[Free,Free]%*%(Hess%*%G - u)[Free]
    # G[Free][G[Free] < 0] <- 0 # projection
    
    ### backtracking
    G_temp <- G
    gradient <- (Hess %*% G - u)[Free]
    Update <- solve(Hess[Free, Free], gradient)
    G_temp[Free] <- G[Free] - alpha*Update
    G_temp[Free][G_temp[Free] < 0] <- 0 # projection
    while (t(G_temp) %*% Hess %*% G_temp - 2 * t(u) %*% G_temp > t(G) %*% Hess %*% G - 2 * t(u) %*% G + 0.5 * t((G_temp-G)[Free]) %*% (2 * gradient)) {
      alpha <- 0.8 * alpha
      print("Shrinking step size")
      G_temp[Free] <- G[Free] - alpha * Update
      G_temp[Free][G_temp[Free] < 0] <- 0 # projection
    }
    if (sum((G_temp-G)^2) < 1e-30){
      cat("Last update has Euclidean norm", sum((G_temp-G)^2))
      cat("Last step size", alpha)
      G <- G_temp
      break
    }
    G <- G_temp
  }
  return(G)
}

## Obtain (g(0), g(T1),..., g(Tn)) from the output G = (g'(0), g'(T1),..., g'(Tn)) 
# of the previous function by computing cumulative sums
# psi = -g is the estimated decreasing score function
neg_dec_score <- function(X, G, theta_pilot = 0) {
  N <- length(X)
  X_sorted <- sort(abs(X - theta_pilot))
  delta_X <- diff(c(0, X_sorted))
  gvalues <- G[-1] + G[-(N+1)]
  gvalues <- 0.5*gvalues*delta_X
  for (j in 2:N) {
    gvalues[j] <- gvalues[j] + gvalues[j-1]
  }
  gvalues <- c(0, gvalues) # gvalues = (g*(0),g*(T1),...,g*(Tn))
  return(gvalues)
}

# Estimated convex loss function G = -phi is obtained by integrating up once again
# G' = g, neg_loglikelihd = (G(0) := 0, G(T1), ..., G(Tn)), no normalization 
neg_loglikelihd <- function(X, G, gvalues, theta_pilot = 0) {
  N <- length(X)
  X_sorted <- sort(abs(X - theta_pilot))
  delta_X <- diff(c(0, X_sorted))
  Negloglikelihd <- gvalues[-(N + 1)] * delta_X + (G[-(N + 1)]/2 + (G[-1] - G[-(N + 1)]) / 6) * delta_X^2
  for (j in 2:N) {
    Negloglikelihd[j] <- Negloglikelihd[j] + Negloglikelihd[j-1]
  }
  Negloglikelihd <- c(0, Negloglikelihd) # neg_loglikelihd = (G(0):=0,G(T1),...,G(Tn))
  return(Negloglikelihd)
}

# True projected score function psi^* for the Cauchy density:
t0 <- newton(function(t) t*sin(t) + cos(t) - 1, 2)$root
x0 <- cot(t0 / 2)
psi_cauchy <- function(x, x0) {
  if (length(x) > 1) return(sapply(x, psi_cauchy, x0 = x0))
  x <- max(min(x, x0), -x0)
  -2*x / (1 + x^2)
}

############# Main script ############
N <- 2000 # sample size
lambda_0 <- 1e-4 # penalization coefficient
noise <- "gaussian"
# noise <- "laplace"
# noise <- "cauchy"
## Generate data
if (noise == "gaussian") {
  X <- rnorm(N, 0, 1)
} else if (noise == "cauchy") {
  X <- rcauchy(N)
} else if (noise == "laplace") {
  X <- rexp(N)
  X <- X * sign(runif(N, -1, 1))
}

# X <- runif(N,-1,1) # U[-1,1]
# X <- rgnorm(N, mu = 0, alpha = 1, beta = 9) # generalized Gaussian
# X <- runif(N, 0.5, 1) + runif(N, 0.5, 1) # triangular density: log-concave; infinite Fisher information
# X <- runif(N, -1, 1) + rnorm(N, 0, 1) # U[-1,1] * N(0,1)

## Add beta noise:
# Z <- rbeta(N, 5, 5) - 0.5
# X <- X + 0.05*Z

theta_pilot <- 0 # oracle location parameter (center of symmetry)
X_sorted <- sort(abs(X - theta_pilot)) # |X_i - theta| values 
delta_X <- diff(c(0, X_sorted))
tailgrid <- (1:(N / 100)) * mean(delta_X)
# For plotting, add additional grid points to the right of Xmax = X_sorted[N]
Xplot <- c(0, X_sorted, tailgrid + X_sorted[N])

## Estimate g = -psi, the negative score function:
G <- symmetric_neg_decr_score_est(X, lambda_0 = lambda_0)
gvalues <- neg_dec_score(X, G)
Negloglikelihd <- neg_loglikelihd(X, G, gvalues)
NStailY <- gvalues[N + 1] + tailgrid * G[N + 1]
NSYplot <- c(gvalues, NStailY)
## Compute G = -phi, the estimated convex loss (negative log-likelihood):
NLLtailY <- Negloglikelihd[N + 1] + tailgrid * gvalues[N + 1] + tailgrid^2 / 2 * G[N + 1]
NLLYplot <- c(Negloglikelihd, NLLtailY)
Xmin <- min(X_sorted)
Xmax <- max(X_sorted)

# True negative score (-psi_0) and negative log-likelihood (-phi_0) functions:
if (noise == "gaussian") {
  NS0Yplot <- Xplot
  NLL0Yplot <- Xplot^2 / 2
} else if (noise == "cauchy") {
  NS0Yplot <- 2 * Xplot/(1 + Xplot^2) # psi_0 is not monotone in this case
  NLL0Yplot <- log(1 + Xplot^2) # offset Negloglikelihd by log(pi) so that Negloglikelihd(0) = 0
} else if (noise == "laplace") {
  NS0Yplot <- sign(Xplot)
  NLL0Yplot <- abs(Xplot) # offset Negloglikelihd by log(2) so that Negloglikelihd(0) = 0
}

## Create plots: 
mydf <- data.frame(X = Xplot, NScore = NSYplot, NLL = NLLYplot, NScore0 = NS0Yplot, NLL0 = NLL0Yplot)

# Negative score function (-psi) plot: NScore (estimated score, blue) and NScore0 (true score, red) and the data (green)
# TODO: In the Cauchy case, also include the projected score psi^*
plotNS.obj <- ggplot(data = mydf) + 
  geom_line(aes(x = X, y = NScore, color = "NShat")) + 
  geom_line(aes(x = X, y = NScore0, color = "NS0")) +
  scale_color_manual(name = NULL, 
                     breaks = c("NShat", "NS0"), 
                     values = c("NShat" = 4, "NS0" = 2), 
                     labels = c(expression(-hat(psi)), expression(-psi[0]))) +
  ylab("negative score function") + xlab(NULL) +
  geom_vline(xintercept = Xmin, linetype = "dotted") + geom_vline(xintercept = Xmax, linetype = "dotted") + # dotted vertical lines marking Xmin and Xmax
  geom_point(data = data.frame(X = X_sorted), aes(x = X, y = 0), color = 3, shape = 4, size = 0.5) +
  theme_minimal() + theme(legend.position = c(0.13, 0.92))
# Original colours: darkblue -> 4; red -> 2; green -> 3

# Negative log-likelihood (-phi) plot:
plotNLL.obj <- ggplot(data = mydf) + geom_line(aes(x = X, y = NLL, color = "NLLhat")) + geom_line(aes(x = X, y = NLL0, color = "NLL0")) +
  scale_color_manual(name = NULL, breaks = c("NLLhat", "NLL0"), values = c("NLLhat" = 4, "NLL0" = 2), labels = c(expression(-hat(phi)), expression(-phi[0]))) +
  ylab("negative log-likelihood") + xlab(NULL) +
  geom_vline(xintercept = Xmin, linetype = "dotted") + geom_vline(xintercept = Xmax, linetype = "dotted") +
  geom_point(data = data.frame(X = X_sorted), aes(x = X, y = 0), color = 3, shape = 4, size = 0.5) +
  theme_minimal() + theme(legend.position = c(0.13, 0.92))

# Arrange two plots side by side:
grid.arrange(plotNS.obj, plotNLL.obj, ncol = 2)

####################################

###### gvalues = (g(0), g(T1),...,g(Tn)))) #######
###### G = (g'(0), g'(T1),...,g'(Tn)) #######

gtheta <- function(theta, X1, X, theta_pilot = 0) {
  X_sorted <- sort(abs(X - theta_pilot))
  knots <- c(0, X_sorted)
  absdelta <- abs(X1 - theta)
  index <- findInterval(absdelta, knots)
  Values <- gvalues[index] # g(left end point)
  distoleft <- absdelta - knots[index]
  Values <- Values + distoleft*G[index] # linear part
  quadpts <- index != N+1
  quadindex <- index[quadpts]
  Values[quadpts] <- Values[quadpts] + 0.5*(G[quadindex+1] - G[quadindex])*distoleft[quadpts]^2/delta_X[quadindex]
  Values <- Values*sign(X1 - theta)
  return(sum(Values))
}

### gprimetheta is dg/dx = -dg/d theta
gprimetheta <- function(theta, X1, X, theta_pilot = 0) {
  X_sorted <- sort(abs(X - theta_pilot))
  knots <- c(0, X_sorted)
  absdelta <- abs(X1 - theta)
  index <- findInterval(absdelta, knots)
  Values <- G[index] # g(left end point)
  quadpts <- index != N+1
  quadindex <- index[quadpts]
  distoleft <- absdelta - knots[index]
  Values[quadpts] <- Values[quadpts] + (G[quadindex+1] - G[quadindex])*distoleft[quadpts]/delta_X[quadindex]
  return(sum(Values))
}

# ####### Compare with the sample mean
# Thetahats <- NULL
# Thetatildes <- NULL
# for (j in 1:1000) {
#   # X1 <- runif(N, -1, 1)
#   # X1 <- rnorm(N)
#   X1 <- rgnorm(N, mu = 0, alpha = 1, beta = 4)
#   thetatilde <- mean(X1)
#   thetahat <- thetatilde # initialization
# 
#   while (abs(gtheta(thetahat, X1, X, theta_pilot)) > 1e-13 * N) {
#     ### backtracking
#     # alpha <- 1
#     # Update <- gtheta(thetahat,X1, X, theta_pilot)/gprimetheta(thetahat, X1, X, theta_pilot)
#     # if(abs(Update) < 1e-12){break}
#     # thetatemp <- thetahat + alpha*Update
#     # while (abs(gtheta(thetatemp,X1, X, theta_pilot)) > abs(gtheta(thetahat, X1, X, theta_pilot))) {
#     #   alpha <- 0.8*alpha
#     #   thetatemp <- thetahat + alpha*Update
#     # }
#     # #print(alpha)
#     # #print(abs(gtheta(thetatemp, X1, X, theta_pilot)) - abs(gtheta(thetahat, X1, X, theta_pilot)))
#     # thetahat <- thetatemp
#     # no backtracking
#     Update <- gtheta(thetahat, X1, X, theta_pilot)/gprimetheta(thetahat, X1, X, theta_pilot)
#     thetahat <- thetahat + Update
#     if(abs(Update) < 1e-12){break}
#   }
#   #print(abs(gtheta(thetahat, X1, X, theta_pilot)))
#   Thetahats <- c(Thetahats, abs(thetahat))
#   Thetatildes <- c(Thetatildes, abs(thetatilde)/sqrt(2))
#   # divide the |mean - theta0| by sqrt(2) to balance the extra 1000 samples used for score matching
# }
# mean(Thetahats)
# mean(Thetatildes)

# prop <- function(x, n){
#   as.numeric(round(prop.test(x = x, n = n, conf.level = 0.95, correct = FALSE)$conf.int * 100, 2))
# }
# 
# binconf <- function(x){
#   cbind(round(x / sum(x) * 100, 2), t(sapply(x, prop, n = sum(x))))
# }
