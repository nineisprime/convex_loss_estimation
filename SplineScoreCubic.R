library(ggplot2)
library(gridExtra)
library(Matrix)

# rm(list = ls())  # clear all variables in workspace
# graphics.off()  # clear all plots

N <- 800 # sample size
lambda0  <- 0.001 # penalization coefficient

# Additive noise
# Noise <- runif(N, -1, 1)
Noise <- rnorm(N, 0, 1)
# Noise <- c(Noise, rnorm(50, 0, 1))
Noise <- rnorm(N, 0, 1)
negNoise <- (Noise < 0)
# Asymmetric log-concave noise variables with mean = median = 0
Noise <- runif(N, -2*sqrt(2/pi),0)*negNoise + Noise*(1-negNoise)

Residualsorted <- sort(Noise)
deltaResidual <- diff(Residualsorted)

lambda <- lambda0*2*N*sd(Noise)^3 # sd(X)^3: to make lambda scale invariant
A <- diag(N)
diag(A[-1,]) <- -1

B <- diag(N)
diag(B[-1,]) <- 1

W <- diag(c(0,1/deltaResidual))

v <- N:1
AAt_inv <- pmin(matrix(v, nrow = N, ncol = length(v), byrow = TRUE), matrix(v, nrow = N, ncol = length(v), byrow = FALSE)) #solve(A%*%t(A))
A_Binv <- matrix(2, ncol = N, nrow = N)
A_Binv[upper.tri(A_Binv, diag = TRUE)] <- 0
A_Binv <- A_Binv + diag(N)
A_Binv <- A_Binv*sign(toeplitz(1:N%%2) - 0.5)#A%*%solve(B)

TEMP <- diag(c(1, deltaResidual^2))%*%AAt_inv + 12*lambda*W
Schur <- TEMP[-1, -1]
Schur <- solve(Schur - matrix(TEMP[-1, 1])%*%t(matrix(TEMP[1, -1]))/TEMP[1, 1])
TEMP <- -t(matrix(TEMP[1, -1]))%*%Schur/TEMP[1, 1]
TEMP <- 6*lambda*rbind(TEMP, Schur)
TEMP <- cbind(rep(0, N), TEMP)
# TEMP is D in the spline note
# gstar = (g*(T1),...,g*(Tn)) is the negative score; gstarprime = (g*'(T1),...,g*'(Tn))
gstarprime <- solve(lambda*(t(A_Binv)%*%W%*%A_Binv + 3*W) - 6*lambda*W^2%*%TEMP, (N:1)%%2)
gstar <- TEMP%*%gstarprime
gstarprime <- cumsum(gstarprime*sign((1:N)%%2 - 0.5))*sign((1:N)%%2 - 0.5)#solve(B)%*%gstarprime
gstar <- cumsum(gstar) #solve(A)%*%gstar
#sum(gstar) # check == 0

############# plot g* ############
# Evaluate a cubic spline at x:
gstars <- function(x,knots,gstar,gstarprime){
  N <- length(knots)
  deltaknots <- diff(knots)
  secantslopes <- diff(gstar)/deltaknots
  sumgprime <- gstarprime[-N] + gstarprime[-1]
  sumgprime1 <- 2*gstarprime[-N] + gstarprime[-1]
  gvalues <- rep(0, length(x))
  index <- findInterval(x, knots)
  gvalues[index == 0] <- gstar[1] + gstarprime[1]*(x[index == 0] - knots[1])
  gvalues[index == N] <- gstar[N] + gstarprime[N]*(x[index == N] - knots[N])
  cubicpts <- !index%in%c(0, N)
  cubicindex <- index[cubicpts]
  y <- (x[cubicpts] - knots[cubicindex])/deltaknots[cubicindex]
  temp <- gstarprime[cubicindex]*y + (3*secantslopes[cubicindex] - sumgprime1[cubicindex])*y^2 + (sumgprime[cubicindex] - 2*secantslopes[cubicindex])*y^3
  temp <- temp*deltaknots[cubicindex] + gstar[cubicindex]
  gvalues[cubicpts] <- temp
  return(gvalues)
}

xgrid <- (-1100:1100)/1000 # grid of x values

tailgrid <- (1:max((N/100),7))*mean(deltaResidual)
tailY <- gstar[N] + tailgrid*gstarprime[N]
Negtailgrid <- (-max((N/100),7):-1)*mean(deltaResidual)
NegtailY <- gstar[1] + Negtailgrid*gstarprime[1]

Xplot <- c(Negtailgrid + Residualsorted[1], Residualsorted, tailgrid + Residualsorted[N], xgrid)
Yplot <- c(NegtailY, gstar, tailY, gstars(xgrid, Residualsorted, gstar, gstarprime))

mydf <- data.frame(X = Xplot, Y = Yplot)
plot0.obj <- ggplot(data = mydf) + geom_line(aes(x = X, y = Y), color = "blue") + geom_line(aes(x = X, y = (sign(X)+1)*X/2), color = "red") +
  geom_point(aes(x = Residualsorted[1], y = gstar[1]), colour = "black") +
  geom_point(aes(x = Residualsorted[N], y = gstar[N]), colour = "black")
grid.arrange(plot0.obj)

# Red line: true score function
# Blue line: score function estimate
