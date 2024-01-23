library(fdrtool)
library(pracma)
library(manipulate)
# library(spatstat.core)

# Compute a decreasing score estimate via antitonic projection
decr_score_est <- function(x, k = 1000, score_est = NULL, kernel = "gaussian", bw = bw.nrd0(x), grid = seq(min(x), max(x), by = IQR(x) / k), kernel_pts = max(min((max(x) - min(x)) / 0.05, 2^10), 2^20), debug = FALSE){
  if (debug) browser()
  # x: data points
  # k: number of points at which to evaluate Jhat and psi^*
  u <- 0:k / k
  
  ## Kernel method
  # Compute a kernel density estimate f, and evaluate its corresponding quantile function Q = F^(-1) together with J = f \circ F^(-1) on u:
  kde <- density(x, kernel = kernel, bw = bw, n = kernel_pts)
  f <- approxfun(kde$x, kde$y, method = "linear", rule = 2)
  Fx <- cumsum(kde$y * c(0, diff(kde$x)))
  Fx <- Fx / Fx[length(kde$x)]
  F <- approxfun(kde$x, Fx, method = "linear", rule = 2)
  # In case F has constant pieces:
  Fx[duplicated(Fx)] <- Fx[duplicated(Fx)] + 1e-12
  Q <- approxfun(Fx, kde$x, method = "linear", rule = 2, ties = max)
  # No need to determine quantiles numerically:
  # Q <- function(u){brent(function(x) {F(x) - u}, min(x) - 1, max(x) + 1)$root}
  # Q <- sapply(u, Q)
  # The quantile.density function from spatstat.core does not always give accurate results:
  # Q <- quantile(kde, probs = u)
  J <- f(Q(u))
  
  # Compute the least concave majorant Jhat of J on u, together with its right derivative:
  lcm <- gcmlcm(u, J, type = "lcm")
  slopes <- lcm$slope.knots
  JhatR <- stepfun(lcm$x.knots, c(slopes[1], slopes, slopes[length(slopes)]))
  # Return the values of psi^* = JhatR \circ F on a grid:
  psi <- JhatR(F(grid))
  phi <- c(0, cumsum(psi[-k] * diff(grid)))
  pdf <- exp(phi); pdf <- pdf / sum((pdf[-1] + pdf[-k]) * diff(grid) / 2); phi <- log(pdf)
  scorefun <- list(obs = x, x = Q, kde = kde, bw = bw, kernel = kernel, u = u, psi = psi, phi = phi, pdf = pdf, J = J, lcm = lcm, grid = grid)
  class(scorefun) <- "scorefun"
  scorefun
}

# psi_0 and psi_0^* for the Cauchy density:
psi_cauchy <- function(x){
  if (length(x) > 1) return(sapply(x, psi_cauchy))
  t0 <- newton(function(t) t * sin(t) + cos(t) - 1, 2)$root
  x0 <- cot(t0 / 2)
  x <- c(x, max(min(x, x0), -x0))
  # Return psi_0(x) and psi_0^*(x)
  -2*x / (1 + x^2)
}

# phi_0 and phi_0^* for the Cauchy density:
phi_cauchy <- function(x){
  if (length(x) > 1) return(sapply(x, phi_cauchy))
  t0 <- newton(function(t) t * sin(t) + cos(t) - 1, 2)$root
  x0 <- cot(t0 / 2)
  # Return l_0(x) and l_0^*(x)
  c(-log(1 + x^2) - log(pi), ifelse(abs(x) <= x0, -log(1 + x^2) - log(pi), -(abs(x) - x0) * sin(t0) - log(1 + x0^2) - log(pi)))
}

# matplot(grid <- seq(-10, 10, length.out = 400), t(psi_cauchy(grid)), type = 'l', lty = c(1, 1), col = c(2, 3), xlab = "x", ylab = expression(psi))
# matplot(grid <- seq(-10, 10, length.out = 400), -t(phi_cauchy(grid)), type = 'l', lty = c(1, 1), col = c(2, 3), xlab = "x", ylab = "y")

# psi <- function(x, ...) UseMethod("psi")
# phi <- function(x, ...) UseMethod("phi")
# 
# psi.scorefun <- function(scorefun, x){
#   stepfun(scorefun$grid, c(scorefun$psi[1], scorefun$psi))(x)
# }

# Plot (true and estimated) decreasing score functions, together with the corresponding loss functions and densities
plot_fisher <- function(scorefun = NULL, dist = FALSE, plot.main = TRUE, k = 1000, debug = FALSE, alpha = 2, sigma = 2, save = TRUE, filename = NULL, ...) {
  if (debug) {
    save <- FALSE; browser()
  }
  kde <- NULL
  empirical <- (class(scorefun) == "scorefun") # whether or not the score function is estimated from data
  if (empirical) list2env(scorefun, envir = environment())
  else stopifnot("scorefun object or named distribution required as input" = !isFALSE(dist))
  if (save & is.null(filename)) filename <- ifelse(is.null(scorefun), dist, paste0(dist, "-data"))
  main <- NULL # plot title
  
  ## J plot:
  if (save) pdf(file = paste0(filename, "-J.pdf"), width = 8.2, height = 5.4)
  if (plot.main) main <- "(a) Plot of J and its least concave majorant"
  # main <- expression(bold(paste("Plot of ", J, " and ", hat(J))))
  legend <- col <- NULL # legend text and colours
  u <- 0:k / k # grid for the J plot
  J0 <- J0hat <- Jhat <- NULL
  # Population-level J_0 and its LCM for some specific densities
  if (dist == "gaussian") J0 <- J0hat <- dnorm(qnorm(u))
  else if (dist == "cauchy") {
    J0 <- (1 - cos(2*pi*u)) / (2*pi)
    lcm0 <- gcmlcm(u, J0, type = "lcm")
    J0hat <- approxfun(lcm0$x.knots, lcm0$y.knots)(u)
    # t0 <- newton(function(t) t * sin(t) + cos(t) - 1, 2)$root; x0 <- cot(t0/2)
    # segments(0, 0, t0 / (2 * pi), (1 - cos(t0)) / (2 * pi), col = 3)
    # segments(1, 0, 1 - t0 / (2 * pi), (1 - cos(t0)) / (2 * pi), col = 3)
  }
  else if (dist == "pareto") {
    stopifnot(alpha > 0 & sigma > 0)
    J0 <- alpha * (2 * pmin(u, 1-u))^(1 + 1/alpha) / (2 * sigma)
    J0hat <- alpha * (2 * pmin(u, 1-u)) / (2 * sigma)
  }
  # LCM based on an estimated J:
  if (empirical) {
    Jhat <- approxfun(lcm$x.knots, lcm$y.knots)(u)
    legend <- c("J", expression(hat(J))); col <- c(1, 4)
    matplot(u, cbind(J, Jhat), col = c(1, 4), lty = c(1, 1), type = "l", ylim = 1.1 * range(c(Jhat, J0hat)), xlab = "u", ylab = "J", main = main, ...)
  } else plot(range(u), 1.1 * range(J0), type = "n", xlab = "u", ylab = "J", main = main, ...)
  if (!is.null(J0)) {
    matlines(u, cbind(J0, J0hat), col = c(2, 3), lty = c(1, 1), lwd = c(2, 2))
    legend <- c(legend, expression(J[0]), expression(hat(J)[0])); col <- c(col, 2, 3)
  }
  legend("topright", legend = legend, col = col, lty = rep(1, 4), ...)
  if (save) dev.off()
  # Default R colours: 1 - black, 2 - red, 3 - green, 4 - blue, 5 - cyan, 6 - magenta, 7 - yellow, 8 - grey
  
  ## Score function (hat(psi)) plot:
  if (save) pdf(file = paste0(filename, "-score.pdf"), width = 8.2, height = 5.4)
  if (plot.main) main <- ("(b) Score function plot")
  grid1 <- psi0 <- psihat <- legend <- col <- NULL
  # Define an appropriate grid (xlim) and legend
  if (empirical) {
    legend <- "Estimated score"; col <- 1
    if (dist == "cauchy") {
      grid1 <- grid[abs(grid) < 3]; psihat <- psi[abs(grid) < 3]
    } else {
      grid1 <- grid; psihat <- psi
    }
  } else {
    if (dist == "pareto") {
      grid1 <- (10^(min(16, 1/alpha)) - 1) * sigma * seq(-1, 1, length.out = k) 
    } else grid1 <- seq(-5, 5, length.out = k)
  }
  # True score functions and their antitonic projections
  if (dist == "gaussian") {
    psi0 <- -grid1; legend <- c(legend, "Gaussian score"); col <- c(col, 2)
  } else if (dist == "cauchy") {
    psi0 <- t(psi_cauchy(grid1)); legend <- c(legend, "Cauchy score", "Projected Cauchy score"); col <- c(col, 2, 3)
  } else if (dist == "pareto") {
    psi0 <- cbind(-(alpha + 1) * sign(grid1) / (abs(grid1) + sigma) , -alpha * sign(grid1) / sigma); 
    legend <- c(legend, paste0("Pareto", "(", alpha, ",", sigma, ") score"), "Projected Pareto score"); col <- c(col, 2, 3)
  }
  # Plot score functions
  plot(range(grid1), 1.1 * range(c(psi0, psihat)), type = "n", xlab = "x", ylab = expression(psi), main = main, ...)
  if (empirical) lines(grid1, psihat, type = "l")
  if (!is.null(psi0)) matlines(grid1, psi0, type = "l", col = c(2, 3), lty = c(1, 1), lwd = c(2, 2))
  legend("topright", legend = legend, col = col, lty = rep(1, 4), ...)
  if (save) dev.off()
  
  ## Loss function (-phi) plot:
  if (save) pdf(file = paste0(filename, "-loss.pdf"), width = 8.2, height = 5.4)
  if (plot.main) main <- ("(c) Loss function plot")
  legend <- col <- l0 <- loss <- NULL
  if (empirical) {
    legend <- c("Estimated convex loss"); col <- 1; loss <- -phi
  } else {
    if (dist == "gaussian") grid <- seq(-5, 5, length.out = k)
    if (dist == "cauchy") grid <- seq(-9, 9, length.out = k)
    if (dist == "pareto") grid <- (10^(min(16, 1/alpha)) - 1) * sigma * seq(-1, 1, length.out = k)
  }
  # Population-level loss functions:
  if (dist == "gaussian") {
    l0 <- -log(dnorm(grid)); legend <- c(legend, "Squared error loss"); col <- c(col, 2)
  } else if (dist == "cauchy") {
    l0 <- -t(phi_cauchy(grid)); legend <- c(legend, "-log(Cauchy density)", "Projected Cauchy loss"); col <- c(col, 2, 3)
  } else if (dist == "pareto") {
    l0 <- -log(cbind(alpha * sigma^alpha / (2 * (abs(grid) + sigma)^(alpha + 1)), # true density
      alpha * exp(-alpha * abs(grid) / sigma) / (2 * sigma))) # Fisher projection
    legend <- c(legend, paste0("-log(Pareto", "(", alpha, ",", sigma, ")", " density)"), "Fisher-projected loss"); col <- c(col, 2, 3)
    if (alpha > 1) {
      l0 <- cbind(l0, -log((alpha - 1) * exp(-(alpha - 1) * abs(grid) / sigma) / (2 * sigma))) # log-concave MLE
      legend <- c(legend, "ML-projected loss"); col <- c(col, 4)
    }
  }
  # Plot loss functions
  plot(range(grid), 1.1 * range(c(l0, loss)), type = "n", xlab = "x", ylab = "y", main = main, ...)
  if (empirical) lines(grid, loss, type = "l")
  if (!is.null(l0)) matlines(grid, l0, type = "l", col = c(2, 3, 4), lty = c(1, 1), lwd = c(2, 2))
  legend("bottomright", legend = legend, col = col, lty = rep(1, 5), ...)
  if (save) dev.off()

  ## Density plot:
  if (save) pdf(file = paste0(filename, "-density.pdf"), width = 8.2, height = 5.4)
  if (plot.main) main <- ("(d) Density plot")
  legend <- col <- p0 <- NULL
  if (!is.null(l0)) p0 <- exp(-l0)
  if (empirical) {
    # legend <- c("Estimated log-concave density", "KDE"); col <- c(1, 6)
    legend <- c("Estimated log-concave density"); col <- 1
  } else pdf <- NULL
  # Use the same grid as for the loss function plot
  # True densities and their projections given by exp(-l0)
  if (dist == "gaussian") {
    legend <- c(legend, "Gaussian density"); col <- c(col, 2)
  } else if (dist == "cauchy") {
    legend <- c(legend, "Cauchy density", "Fisher projection"); col <- c(col, 2, 3)
  } else if (dist == "pareto") {
    legend <- c(legend, paste0("Pareto", "(", alpha, ",", sigma, ")", " density"), "Fisher projection"); col <- c(col, 2, 3)
    if (alpha > 1) {
      legend <- c(legend, "Log-concave ML projection"); col <- c(col, 4)
    }
  }
  # Plot densities
  plot(range(grid), 1.1 * range(c(p0, pdf)), type = "n", xlab = "x", ylab = "y", yaxs = "i", main = main, ...)
  if (empirical) {
    lines(grid, pdf)
    # lines(kde, col = 6)
  } 
  if (!is.null(l0)) matlines(grid, p0, type = "l", col = c(2, 3, 4), lty = rep(1, 3), lwd = c(2, 2))
  legend("topright", legend = legend, col = col, lty = rep(1, 6), ...)
  if (save) dev.off()
  # ## KDE plot (based on data):
  # if (empirical) {
  #   if (save) pdf(file = paste0(filename, "-kde.pdf"), width = 8.2, height = 5.4)
  #   if (plot.main) main <- "Kernel density estimate"
  #   plot(kde, xlim = c(min(grid), max(grid)), main = main)
  #   if (save) dev.off()
  # }
}

d_gmixture <- function(x, mu, p = 0.5) (1 - p) * dnorm(x + mu) + p * dnorm(x - mu)
p_gmixture <- function(x, mu, p = 0.5) (1 - p) * pnorm(x + mu) + p * pnorm(x - mu)

quantile_gmixture <- function(u, mu, p = 0.5){
  if (length(u) > 1) return(sapply(u, quantile_gmixture, mu = mu, p = p))
  brent(function(x){(1 - p) * pnorm(x + mu) + p * pnorm(x - mu) - u}, qnorm(u, -mu, 1), qnorm(u, mu, 1), tol = 1e-12)$root
}

J_gmixture <- function(u, mu, p = 0.5){
  if (length(u) > 1) return(sapply(u, J_gmixture, mu = mu, p = p))
  d_gmixture(quantile_gmixture(u, mu, p), mu, p)
}

set.seed(2)

## Pareto example
plot_fisher(dist = "pareto", alpha = 2, sigma = 3)

n <- 2000
## Cauchy example
x <- rcauchy(n)
# manipulate(plot(density(x, bw = bw, kernel = kernel), 
#                 main = paste0("Optimal bandwidths: h_NRD = ", round(bw.nrd0(x), 3), "; h_CV = ", round(bw.ucv(x), 3))), 
#            bw = slider(1/n, 1, step = 1/n), kernel = picker("gaussian", "epanechnikov", "triangular"))
score_est <- decr_score_est(x, kernel_pts = 2^16, grid = seq(-10, 10, length.out = 1000))
plot_fisher(scorefun = score_est, dist = "cauchy")
# Population-level plots
plot_fisher(dist = "cauchy")

n <- 1e+5
## Gaussian example
x <- rnorm(n)
score_est <- decr_score_est(x, kernel_pts = 2^16)
plot_fisher(score_est, dist = "gaussian", cex = 0.9)

## Gaussian example with an undersmoothed KDE (small bandwidth)
x <- rnorm(n)
score_est <- decr_score_est(x, kernel_pts = 2^16, bw = 0.01)
plot_fisher(score_est, dist = "gaussian", cex = 0.9, save = FALSE)

## Asymmetric example
p <- 1/(1 + sqrt(2/pi))
epsilon <- rbinom(n, 1, p)
x <- -abs(rnorm(n)) * (1 - epsilon) + abs(rcauchy(n)) * epsilon
grid <- seq(-3, 10, length.out = 1000)
score_est <- decr_score_est(x, kernel_pts = 2^16, grid = grid)
plot_fisher(score_est, filename = "asymmetric", save = FALSE)
lines(grid[grid < 0], 2 * (1-p) * dnorm(grid[grid < 0]), col = 2)
lines(grid[grid >= 0], 2 * p * dcauchy(grid[grid >= 0]), col = 2)
# TODO: Add red lines to the J, score and loss function plots

## Mixture of uniforms: 3/4 * U[-1,0] + 1/4 * U[0,1]
epsilon <- rbinom(n, 1, 1/4)
x <- -runif(n) * (1 - epsilon) + runif(n) * epsilon
score_est <- decr_score_est(x, kernel_pts = 2^16)
plot_fisher(score_est, filename = "uniform-mixture")

proj_cauchy <- function(x){
  if (length(x) > 1) return(sapply(x, proj_cauchy))
  t0 <- newton(function(t) t * sin(t) + cos(t) - 1, 2)$root
  x0 <- cot(t0/2)
  # x0 <- brent(function(x) 2 * x * atan(1/x) - 1, 0, 1)$root
  if (abs(x) <= x0) return(1 / (pi * (1 + x^2)))
  else exp(-sin(t0) * (abs(x) - x0)) / (pi * (1 + x0^2))
}

huber <- function(x, K){
  if (length(x) > 1) return(sapply(x, huber, K))
  ifelse(abs(x) < K, x^2/2, K*abs(x) - K^2/2)
}

# Plot of the Huber and squared error loss
plot(x <- seq(-5, 5, by = 0.01), huber(x, K = 2), type = 'l', col = '2', lwd = 2, xlim = c(-5, 5), ylim = c(-0.2, 25/2),
     xlab = '', ylab = '', axes = FALSE, xaxs = 'i', yaxs = 'i')
axis(1, at = c(-6, -2, 2, 6), pos = 0, labels = c("", "-K", "K", ""))
segments(0, 0, 0, 25/2)
# abline(h = 0)
abline(v = -2, lty = 3)
abline(v = 2, lty = 3)
# axis(2, labels = F, pos = 0)
arrows(5 - 0.1, 0, 5, 0, angle = 30, length = 0.1)  # Arrow for x-axis
arrows(0, 25/2 - 0.1, 0, 25/2, angle = 30, length = 0.1)  # Arrow for y-axis
# axis(2)
lines(x, x^2/2, type = 'l', col = '4')
legend(x = -5.2, y = 2, legend = c("Huber loss", "Squared error loss"), cex = 0.9, lty = c(1, 1), col = c(2, 4), bty = 'n')

# Cauchy log-likelilhood
n <- 25; scale <- 0.01
set.seed(47)
# set.seed(3)
x <- rcauchy(n)
# manipulate(plot(y <- seq(-100, 100, length.out = 1e+4), sapply(y, function(z){-log(pi * scale) - mean(log(1 + ((z - scale * x) / scale)^2))}), 
#      type = 'l', xlab = '', ylab = ''), scale = slider(0.5, 2, step = 0.1))
plot(y <- seq(-10, 10, length.out = 1e+4), sapply(y, function(z){-log(pi * scale) - mean(log(1 + ((z - x) / scale)^2))}), 
     type = 'l', xlab = '', ylab = '', main = 'Plot of the Cauchy log-likelihood')
# plot(y <- seq(-100, 100, length.out = 1e+4), sapply(y, function(z){-log(pi * scale) - mean(log(1 + ((z - x) / scale)^2))}), 
     # type = 'l', xlab = '', ylab = '', main = 'Plot of the Cauchy log-likelihood')
points(x, rep(-10.5, length(x)), pch = 4)
# points(x, rep(6, length(x)), pch = 4)
