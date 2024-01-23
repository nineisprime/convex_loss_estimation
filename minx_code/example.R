
source("scoreFun.R")

n = 500

test1 = FALSE
test2 = TRUE
test3 = FALSE

if (test1){

    Z = rcauchy(n)
    ##Z = rnorm(n)
    ##Z = runif(n, -1, 1)

    out = antitonicScoreMatch(Z, lambda=.5)

    Z_sorted = sort(Z)

    plotLossScore(psi_evals=out$psistars, pts_evals=Z_sorted, xlim=c(-5,5),
                  psi_true_evals= -2*Z_sorted/(1+Z_sorted^2),
                  plot.loss=FALSE)
}

if (test2){

    ##Z = runif(n, -1, 1)
    Z = rcauchy(n)
    kde <- density(Z, n=2^15)
    print(kde$bw)
    out = isotonizeDensity(kde$y, kde$x)

    Z_sorted = sort(Z)

    plotLossScore(psi_evals=out$psi(Z_sorted), pts_evals=Z_sorted,
                  xlim=c(-5,5),
                  psi_true_evals=-2*Z_sorted/(1+Z_sorted^2), plot.loss=FALSE)
}


if (test3) {

    ##Z = rcauchy(n)
    Z = rnorm(n)
    out = splineScoreEst(Z, lambda=0.5)

    Z_sorted = sort(Z)

    out2 = isoProjInterpolate(out$psistars, Z_sorted)

    plotLossScore(psi_evals=out2$iso_fn(Z_sorted), pts_evals=Z_sorted,
                  xlim=c(-3,3),
                  psi_true_evals=-2*Z_sorted/(1+Z_sorted^2), plot.loss=FALSE)
}
