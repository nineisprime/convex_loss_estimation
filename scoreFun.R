getPhiPDF <- function(psi_evals, pts_evals) {
    K = length(pts_evals)
    phi <- c(0, cumsum((psi_evals[-1] + psi_evals[-K]) / 2 * diff(pts_evals)))

    pdf <- exp(phi)
    pdf <- pdf / sum((pdf[-1] + pdf[-K]) * diff(pts_evals) / 2)
    phi <- log(pdf)
}

plotLossScore <- function(psi_evals, pts_evals, psi_true_evals=NULL,
                          loss_true_evals=NULL, xlim=NULL, ylim.psi=NULL,
                          ylim.loss=NULL,
                          length.grid=500, plot.loss=TRUE){
  if (is.null(xlim)){
    xlim = c(min(pts_evals), max(pts_evals))
  }
  
  grid = seq(xlim[1], xlim[2], length.out=length.grid)
  psi_grid = approxfun(pts_evals, psi_evals, rule=2)(grid)
  
  if (!is.null(psi_true_evals)){
    psi_true_grid = approxfun(pts_evals, psi_true_evals, rule=2)(grid)
    if (is.null(loss_true_evals))
      loss_true_grid = -getPhiPDF(psi_true_grid, grid)
  }
  
  if (is.null(ylim.psi)){
    delta = 0.1 * diff(range(c(psi_grid, psi_true_grid)))
    ylim.psi = range(c(psi_grid, psi_true_grid)) + c(-delta, delta)
  }
  
  matplot(grid, cbind(psi_grid, psi_true_grid), col=c(1,4), lty=c(1,1),
          type="l", ylim=ylim.psi, xlim=xlim)
  
  if (plot.loss){
    if (!is.null(loss_true_evals)){
      loss_true_grid = approxfun(pts_evals, loss_true_evals, rule=2)(grid)
    }
    
    loss_grid = -getPhiPDF(psi_grid, grid)
    
    if (is.null(ylim.loss)){
      delta = 0.1 * diff(range(c(loss_grid, loss_true_grid)))
      ylim.loss = range(c(loss_grid, loss_true_grid)) + c(-delta, delta)
    }
    
    matplot(grid, cbind(loss_grid, loss_true_grid), col=c(1,4), lty=c(1,1),
            type="l", ylim=ylim.loss, xlim=xlim)
  }
}

## Gives
##   iso_fn = IsoRegress(init_fn)
## represented as a function, and its derivative
##
## Assumes "init_fn" are evaluations on "grids"
##
##
## INPUT:  init_fn -- vector of evaluations
##         grids -- points in R at which we evaluate the initial score
##
## OUTPUT: iso_fn -- function  R -> R
##         iso_fn_deriv -- function R -> R
## note that if grids is not the quantiles, the output is not the antitonic projection
isoProjInterpolate <- function(init_fn, grids, return_fcn = TRUE){
  library(Iso)
  library(stats)
  
  K = length(grids)
  
  iso_fn_evals <- pava(init_fn, decreasing = TRUE,
                       long.out = FALSE, stepfun = FALSE)

  non_redundent_index <- diff(iso_fn_evals) != 0
  temp <- c(non_redundent_index, FALSE)
  non_redundent_index <- c(TRUE, non_redundent_index)
  non_redundent_index <- non_redundent_index | temp
  
  
  iso_fn_deriv_evals <- diff(iso_fn_evals)/diff(grids)
  grids1 <- (grids[-K] + grids[-1]) / 2
  
  non_redundent_index1 <- diff(iso_fn_deriv_evals) != 0
  temp1 <- c(non_redundent_index1, FALSE)
  non_redundent_index1 <- c(TRUE, non_redundent_index1)
  non_redundent_index1 <- non_redundent_index1 | temp1
  
  if(return_fcn){
    iso_fn <- approxfun(grids[non_redundent_index], iso_fn_evals[non_redundent_index], method = "linear",rule=2)
    # use linear interpolation for the derivative instead of piece-wise constant
    iso_fn_deriv <- approxfun(grids1[non_redundent_index1], iso_fn_deriv_evals[non_redundent_index1], method = "linear",rule=2)
    return(list(iso_fn=iso_fn, iso_fn_deriv=iso_fn_deriv))
  }else {
    return(list(iso_fn_evals=iso_fn_evals[non_redundent_index], fn_eval_grid=grids[non_redundent_index],
                iso_fn_deriv_evals=iso_fn_deriv_evals[non_redundent_index1], fn_deriv_eval_grid=grids1[non_redundent_index1]))
  }

  
  # return(list(iso_fn = iso_fn_evals, iso_fn_x = grids,
  #             iso_fn_deriv = diff(iso_fn_evals)/diff(grids), iso_fn_deriv_x = (grids[-K] + grids[-1]) / 2))
}

## Gives
## LCM( densities cdot Q )^{(R)} cdot F
## where Q, F are the quantile and CDF with respect to
## "densities", represented as evaluations over "grids"
##
##
## INPUT: densities -- vector of density values
##        grids -- domain of "densities"
##
##        k -- number of discrete points between [0,1] at which
##             to evaluate Jhat and (psi^*)??
##
## OUTPUT: psi -- function  R -> R
##         psi_deriv -- function  R -> R
##
isotonize_score_given_density <- function(grids, densities, k=1000){
  library(fdrtool)

  K <- length(grids)
  
  # Interpolate between density values on the grid
  f <- approxfun(grids, densities, method = "linear", rule = 2)
  Fx <- c(0, cumsum((densities[-1] + densities[-K])/2 * diff(grids)))
  Fx <- Fx / Fx[K]
  # # In case F has constant pieces:
  # Fx[duplicated(Fx)] <- Fx[duplicated(Fx)] + 1e-12
  Q <- approxfun(Fx, grids, method = "linear", rule = 2, ties = max)
  u <- 0:k / k
  J <- f(Q(u))
  init_scores <- diff(J) * k #psihat\circle Qhat
  res = isoProjInterpolate(init_scores, (Q(u[-1]) + Q(u[-(k+1)])) / 2, return_fcn = FALSE)
  psi = approxfun(res$fn_eval_grid, res$iso_fn_evals, method = "linear",rule=2)
  # use linear interpolation for the derivative instead of piece-wise constant
  psi_deriv = approxfun(res$fn_deriv_eval_grid, res$iso_fn_deriv_evals, method = "linear",rule=2)
  
  return(list(psi=psi, psi_deriv=psi_deriv))
}

score_given_density <- function(grids, densities, k=1000){
  library(fdrtool)
  K <- length(grids)
  # Interpolate between density values on the grid
  f <- approxfun(grids, densities, method = "linear", rule = 2)
  Fx <- c(0, cumsum((densities[-1] + densities[-K])/2 * diff(grids)))
  Fx <- Fx / Fx[K]
  Q <- approxfun(Fx, grids, method = "linear", rule = 2, ties = max)
  u <- 0:k / k
  J <- f(Q(u))
  init_scores <- diff(J) * k #psihat\circle Qhat
  psi = approxfun((Q(u[-1]) + Q(u[-(k+1)])) / 2, init_scores, method = "linear",rule=2)
  
  return(psi)
}

############# Method A: kde based #############
#### truncation not implemented
kde_decr_score_est <- function(residuals, k=1000, kernel="gaussian", kernel_pts=2^21,
                               truncation_lim=NULL, set_to_zero=FALSE){
  library(fdrtool)
  library(Iso)
  library(stats)
  residuals <- sort(residuals)
  kde <- density(residuals, kernel = kernel, n = kernel_pts)
  # cat("The largest abs score = ", max(abs(init_scores)), '\n')
  # ####NOT exactly the truncation described in the paper
  # if(!is.null(truncation_lim)){
  #   if(!set_to_zero){
  #     init_scores <- pmax(pmin(init_scores, truncation_lim), -truncation_lim)
  #   }else{
  #     init_scores[(init_scores > truncation_lim) | (init_scores < -truncation_lim)] <- 0
  #   }
  # }
  return(isotonize_score_given_density(kde$x, kde$y, k=max(k, 2*length(residuals))))
}


############# Method B: spline based #############
## Gives
## general spline score matching function and its derivative evaluated at the knots
## no monotone constraint, no symmetricity of the noise assumed
## INPUT:  Z -- vector of the noise proxies, i.e. the residuals
##         lambda0 -- penalization coefficient
##
## OUTPUT: psistars -- (psi*(T1),...,psi*(Tn))
##         psistarprimes -- (psi*'(T1),...,psi*'(Tn))
## where Ti = Z_{(i)}
##
splineScoreEst <- function(Z, lambda0) {
  N <- length(Z)
  lambda <- N * lambda0
  # lambda <- lambda0 * N * sd(Z)^3 #sd(Z)^3: to make lambda scale invariant
  Ts <- sort(Z)
  deltas <- diff(Ts)

  A <- diag(N)
  diag(A[-1,]) <- -1

  B <- diag(N)
  diag(B[-1,]) <- 1

  W <- diag(c(0,1/deltas))

  v <- N:1
  AAt_inv <- pmin(matrix(v, nrow = N, ncol = length(v), byrow = TRUE), 
                  matrix(v, nrow = N, ncol = length(v), byrow = FALSE)) #solve(A%*%t(A))
  A_Binv <- matrix(2, ncol = N, nrow = N)
  A_Binv[upper.tri(A_Binv, diag = TRUE)] <- 0
  A_Binv <- A_Binv + diag(N)
  A_Binv <- A_Binv*sign(toeplitz(1:N%%2) - 0.5)#A%*%solve(B)

  TEMP <- diag(c(1, deltas^2))%*%AAt_inv + 12*lambda*W
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

  return(list("psistars" = -gstar, "psistarprimes" = -gstarprime))
}
## Gives
## evaluations of the general spline score matching function on xs
##
## INPUT:  xs -- evaluation points
##         knots -- (Ti)
##         psistars -- (psi*(T1),...,psi*(Tn))
##         psistarprimes -- (psi*'(T1),...,psi*'(Tn))
##
## OUTPUT: psixs -- psi*(xs)
##         
##
splineScore <- function(xs, knots, psistars, psistarprimes){
  N<-length(knots)
  deltaknots<-diff(knots)
  secantslopes<-diff(psistars)/deltaknots
  sumpsiprime<-psistarprimes[-N]+psistarprimes[-1]
  sumpsiprime1<-2*psistarprimes[-N]+psistarprimes[-1]
  psixs<-rep(0, length(xs))
  index<-findInterval(xs, knots)

  int_indicator = index > 0 & index < N
  int_xs = xs[int_indicator]
  int_index = index[int_indicator]

  psixs[index==0] <- psistars[1] + psistarprimes[1] * (xs[index==0] - knots[1])
  psixs[index==N] <- psistars[N] + psistarprimes[N] * (xs[index==N] - knots[N])

  y <- (int_xs - knots[int_index]) / deltaknots[int_index]
  temp <- psistarprimes[int_index] * y + (3 * secantslopes[int_index] - sumgprime1[int_index]) * y^2 +
    (sumgprime[int_index] - 2 * secantslopes[int_index]) * y^3
  temp<-temp*deltaknots[int_index]+psistars[int_index]
  psixs[int_indicator]<-temp
  return(psixs)
}


## monotone spline without the symmetricity assumption
## INPUT:  Z -- vector of the noise proxies, i.e. the residuals
##         lambda -- penalization coefficient
##
## OUTPUT: psistars -- (psi*(T1),...,psi*(Tn))
##         psistarprimes -- (psi*'(T1),...,psi*'(Tn))
## We only implement the version that Ti=Z_{(i)}.
## The argimin is piecewise quadratic and is uniquely determined by psistars and psistarprimes.
## The evaluation of psistar can be achieved via eval_quad_Score
antitonicScoreMatch <- function(Z, lambda0, issorted=FALSE){
    library(MASS)
    N <- length(Z)
    lambda <- N * lambda0
    # lambda <- lambda0 * N * sd(Z)^3 #sd(Z)^3: to make lambda scale equivariant
    if(!issorted){
      Z <- sort(Z)
    }
    deltas <- diff(Z)

    u = c(0, rep(1, N))
    L = t(matrix(replicate(N, c(2, deltas)), ncol=N))
    L[upper.tri(L)] = 0
    A<-diag(x = 1, nrow = N, ncol = N+1)
    diag(A[-1, 3:(N+1)]) = 1
    A[1, 1] = -1

    D = diag(1/deltas*lambda)
    B = diag(x = 1, nrow = N-1, ncol = N)
    diag(B[, -1]) = -1
    #B above is not the B in the spline_opt.tex; we removed the first all 0 column here; see the tex file;
    #we only use it to calculate the non-zero lower right block of the second term of the Hessian matrix

    Q = t(B) %*% D %*% B  #non-zero lower right block of the second term of the Hessian matrix
    Q = cbind(rep(0, N), Q)
    Q = rbind(rep(0, N+1), Q) #the second term of the Hessian matrix
    Q = 0.25 * crossprod(L %*% A) + Q #Hessian

    #G=(-g(T_1),g'(T_1),g'(T_2),...,g'(T_N))
    G = solve(Q, u) #initialization
    # in theory, Q is invertible
    # G <- tryCatch({
    #   solve(Q, u)
    # }, error = function(err) {
    #   solve(Q + diag(nrow(Q)), u)
    # })

    #gvals_before = L %*% A %*% G /2 #the unconstrained piecewise quadratic argmin

    G[G<0] <- 0 #projection

    MAXITER = 20
    for (i in 1:MAXITER){
        # cat("Newton loop loss",
        #     t(G)%*%Q%*%G - 2*t(u) %*% G,"\n")
        alpha = 1

        Free = (G > 0) | (Q %*% G - u <= 0)

        Gtmp = G
        gradient = (Q %*% G - u)[Free]
        # update = solve(Q[Free, Free], gradient)
        # update <- tryCatch({
        #   solve(Q[Free, Free], gradient)
        # }, error = function(err) {
        #   # ginv(Q[Free, Free])%*%gradient
        #   solve(Q[Free, Free] + diag(nrow(Q[Free, Free])), gradient)
        # })
        update <- solve(Q[Free, Free] + diag(nrow(Q[Free, Free])), gradient)
        Gtmp[Free] = G[Free] - alpha * update
        Gtmp[Free][Gtmp[Free]<0]<-0 #projection

        while (t(Gtmp) %*% Q %*% Gtmp - 2 * t(u) %*% Gtmp >
                t(G) %*% Q %*% G - 2 * t(u) %*% G + 0.5 * t((Gtmp - G)[Free]) %*% (2 * gradient)) {
            alpha <- 0.8 * alpha
            Gtmp[Free] <- G[Free] - alpha * update
            Gtmp[Free][Gtmp[Free]<0] <- 0 #projection
        }
        if (sum((Gtmp - G)^2) < 1e-20 * N){
            # cat("Spline based psihat: Last update's Euclidean norm ", sum((Gtmp-G)^2)^0.5, '\n')
            G = Gtmp
            break
        }
        G = Gtmp
        if(i==MAXITER) {
          cat("Spline based psihat: Last update's Euclidean norm ", sum((Gtmp - G)^2)^0.5, '\n')
        }
    }

    gvals_after = L %*% A %*% G / 2 #gvals_after = (g*(T1),...,g*(Tn))
    # gvalues<-G[-1]
    # gvalues<-gvalues[-1]+gvalues[-N]
    # gvalues<-0.5*gvalues*deltas
    # for (j in 2:(N-1)) {
    #   gvalues[j]<-gvalues[j]+gvalues[j-1]
    # }
    # gvalues<-c(-G[1],-G[1]+gvalues) #gvalues = (g*(T1),...,g*(Tn))
    # gprimes<-G[-1] #gprimes = (g'(T_1),g'(T_2),...,g'(T_N))

    return(list("psistars" = -gvals_after, "psistarprimes" = -G[-1]))
}


## score (and optionally its derivative) evaluated at a set of points, xs
## monotone piecewise quadratic spline without the symmetricity assumption
##
## INPUT:  xs -- points of inquiry
##         psivals -- psi(knots)
##         psiprimes -- psiprime(knots)
##         deltas -- diff(knots)
##
## OUTPUT: psixs -- psi(xs)
##         psiprimexs -- psiprime(xs), if derivative = TRUE
eval_quad_Score <- function(xs, knots, psivals, psiprimes, 
                            robust = FALSE, deltaknots = NULL, derivative = FALSE){

    N = length(knots)
    if (is.null(deltaknots)){
      deltaknots = diff(knots)
    }
    psixs = rep(0, length(xs))

    index = findInterval(xs, knots)

    ## xs inside the range of knots
    int_indicator = index > 0 & index < N
    int_xs = xs[int_indicator]
    int_index = index[int_indicator]

    psivals_left = psivals[int_index]
    distoleft = int_xs - knots[int_index]

    int_psixs = psivals_left + distoleft * psiprimes[int_index] +
        0.5*(psiprimes[int_index+1]-psiprimes[int_index])*distoleft^2 / deltaknots[int_index]

    psixs[int_indicator] = int_psixs

    left_indicator = index == 0
    right_indicator = index == N
    if(!robust){
      psixs[left_indicator] = psivals[1] + (xs[left_indicator] - knots[1]) * psiprimes[1]
      psixs[right_indicator] = psivals[N] + (xs[right_indicator] - knots[N]) * psiprimes[N]
    } else {
      ## set psi(x) for x outside the knots: robust psi
      psixs[left_indicator] = psivals[1]
      psixs[right_indicator] = psivals[N]
    }

    if(derivative){
      if(!robust){
        psiprimexs <- approx(knots, psiprimes, xout = xs, method = "linear", rule = 2)$y
      } else {
        ## set psi'(x) for x outside the knots: robust psi
        psiprimexs <- approx(knots, psiprimes, xout = xs, method = "linear", yleft = 0, yright = 0)$y
      }
      # psiprimexs = rep(0, length(xs))
      # psiprimevals_left = psiprimes[int_index]
      # 
      # int_psiprimexs = psiprimevals_left +
      #   (psiprimes[int_index+1]-psiprimes[int_index])*distoleft / deltaknots[int_index]
      # 
      # psiprimexs[int_indicator] = int_psiprimexs
      # if(!robust){
      #   psiprimexs[left_indicator] = psiprimes[1]
      #   psiprimexs[right_indicator] = psiprimes[N]
      # } else {
      #   ## set psi'(x) for x outside the knots: robust psi
      #   psiprimexs[left_indicator] = 0
      #   psiprimexs[right_indicator] = 0
      # }
      return(list("psixs" = psixs, "psiprimexs" = psiprimexs))
    }
    return(psixs)
}


############# Method B: spline SCALE #############
## evaluate the gradient and the hessian for the linear regression case
psi_psiprime_beta <- function(betahat, X, Y, knots, psivals, psiprimes, 
                              symmetric=FALSE, robust = FALSE, deltaknots = NULL, derivative = TRUE){
  #deltaknots == diff(knots)
  X <- cbind(rep(1, length(Y)), X)
  residuals <- Y - X %*% betahat
  if (derivative) {
    temp <- eval_quad_Score(residuals, knots, psivals, psiprimes, 
                            robust = robust, deltaknots = deltaknots, derivative = derivative)
    psibeta <- t(X) %*% temp[["psixs"]]
    psiprimebeta <- temp[["psiprimexs"]]
    psiprimebeta <- -t(X) %*% diag(psiprimebeta) %*% X
    return(list("psibeta" = psibeta, "psiprimebeta" = psiprimebeta))
  } else {
    temp <- eval_quad_Score(residuals, knots, psivals, psiprimes, 
                            robust = robust, deltaknots = deltaknots, derivative = derivative)
    return(t(X) %*% temp)
  }
}
## Using Newton's method to solve the estimation equations to get betahat
spline_linear_regression_newton <- function(betainit, X, Y, knots, psivals, psiprimes, n_iter = 40,
                                            robust = FALSE, deltaknots = NULL){
  library(MASS)
  if (is.null(deltaknots)){
    deltaknots = diff(knots)
  }
  N <- length(Y)
  betahat <- betainit
  temp <- psi_psiprime_beta(betahat, X, Y, knots, psivals, 
                            psiprimes, robust = robust, deltaknots = deltaknots)
  H <- temp[["psiprimebeta"]]
  v <- temp[["psibeta"]]
  for (l in 1:n_iter) {
    alpha <- 1
    if(sum(v^2)^0.5 < (d + 1) * 1e-6 / N^0.5) {
      # cat("Final l2 norm of the gradient", sum(v^2)^0.5, "\n")
      return(betahat)
    }

    # Update <- ginv(H)%*%v
    Update <- solve(H + diag(nrow(H)), v)
    
    betatemp <- betahat - alpha * Update
    # line search
    # temp <- psi_psiprime_beta(betatemp, X, Y, knots, psivals, 
    #                           psiprimes, robust = robust, deltaknots = deltaknots)
    # #Htemp <- temp[["psiprimebeta"]]
    # vtemp <- temp[["psibeta"]]
    vtemp <- psi_psiprime_beta(betatemp, X, Y, knots, psivals, 
                              psiprimes, robust = robust,
                              deltaknots = deltaknots, derivative = FALSE)
    while (sum(vtemp^2) > 1 * sum(v^2)){
      # print("shrink the step size")
      alpha <- 0.8 * alpha
      betatemp <- betahat - alpha * Update
      # temp <- psi_psiprime_beta(betatemp, X, Y, knots, psivals, 
      #                            psiprimes, robust = robust, deltaknots = deltaknots)
      # vtemp <- temp[["psibeta"]]
      vtemp <- psi_psiprime_beta(betatemp, X, Y, knots, psivals, 
                                 psiprimes, robust = robust,
                                 deltaknots = deltaknots, derivative = FALSE)
    }
    betahat <- betatemp
    temp <- psi_psiprime_beta(betatemp, X, Y, knots, psivals, 
                              psiprimes, robust = robust, deltaknots = deltaknots)
    H <- temp[["psiprimebeta"]]
    v <- vtemp
    print("iter counter")
    print(l)
  }
  cat("Reached the max iter; final l2 norm of the Z function value", sum(vtemp^2)^0.5, "\n")

  return(betahat)
}

###root finding####



# L^1 loss
mae_loss <- function(beta, X, Y) {
  predictions <- X %*% beta
  mean_abs_error <- mean(abs(predictions - Y))
  return(mean_abs_error)
}

max_loss <- function(beta, X, Y) {
  predictions <- X %*% beta
  loss <- max(abs(predictions - Y))
  return(loss)
}