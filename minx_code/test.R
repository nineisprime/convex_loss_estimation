
## Simple experimental setting
N = 400
d = 1

X = matrix(rnorm(N*d), N, d)

beta0 = runif(d) + 1
mu = 1

mask = runif(N) > 0.5
Z = mask*rnorm(N, -2, 0.7) + (1-mask)*rnorm(N, 2, 0.7)
Z = 2*runif(N) - 1

Y = X %*% beta0 + mu + Z

lambda = 0.001 * N

## compute pilot estimator

Xtmp = cbind(rep(1, N), X)
beta_pilot = solve(crossprod(Xtmp), crossprod(Xtmp, Y))

Ztilde_sorted = sort(Y - Xtmp %*% beta_pilot)
deltas = diff(Ztilde_sorted)

u = c(0, rep(1, N))
L = t(matrix(replicate(N, c(2, deltas)), ncol=N))
L[upper.tri(L)] = 0
A<-diag(x = 1, nrow = N, ncol = N+1)
diag(A[-1, 3:(N+1)]) = 1
A[1, 1] = 1


D = diag(1/deltas*lambda)
B = diag(x = 1, nrow = N-1, ncol = N)
diag(B[, -1]) = -1
#D and B above are not the B and D in the tex file;
#we only use them to calculate the non-zero lower right block of the Hessian matrix

Q = t(B) %*% D %*% B  #non-zero lower right block of the Hessian matrix
Q = cbind(rep(0, N), Q)
Q = rbind(rep(0, N+1), Q)
Q = 0.25 * crossprod(L %*% A) + Q

G = solve(Q, u)

gvals_before = L %*% A %*% G

projShape = function(w){
    w[-1][ w[-1] < 0 ] = 0
    return(w)
}

G = projShape(G)

MAXITER = 20
for (i in 1:MAXITER){
    cat("Newton loop loss", 
        t(G)%*%Q%*%G - 2*t(u) %*% G,"\n")
    alpha = 1

    Free = (G > 0) | (Q %*% G - u <= 0)
    Free[1] = TRUE

    Gtmp = G 
    gradient = (Q %*% G - u)[Free]
    update = solve(Q[Free, Free], gradient)

    Gtmp[Free] = G[Free] - alpha * update
    Gtmp[Free] = projShape(Gtmp[Free])

    while (t(Gtmp)%*% Q %*%Gtmp-2*t(u)%*%Gtmp > 
            t(G)%*% Q %*%G-2*t(u)%*%G+0.5*t((Gtmp-G)[Free])%*%(2*gradient)) {
        alpha = 0.8 * alpha
        Gtmp[Free] = G[Free] - alpha*update
        Gtmp[Free] = projShape(Gtmp[Free])
    }
    if (sum((Gtmp-G)^2)<1e-30){
        cat("Last update has Euclidean norm ", sum((Gtmp-G)^2), '\n')
        G = Gtmp
        break
    }
    G = Gtmp
}

gvals_after = L %*% A %*% G



xs = seq(-6, 6, by = 0.01)

max_z = max(Ztilde_sorted)
min_z = min(Ztilde_sorted)
rg_z = max_z - min_z

left_bd = min_z - rg_z*0.1
right_bd = max_z + rg_z*0.1

knots = c(left_bd, Ztilde_sorted, right_bd)

knot_gprimes = c(0, G[-1], 0)

gvals_l = gvals_after[1] - knot_gprimes[2]*(min_z - left_bd)/2 
gvals_r = (knot_gprimes[N+1])*(right_bd - max_z)/2 + gvals_after[N]
knot_gvals = c(gvals_l, gvals_after, gvals_r)

res = evalScore(xs, knots, knot_gvals, knot_gprimes)
# print(res)
res2 = evalScoreDeriv(xs, knots, knot_gvals, knot_gprimes)
# plot(xs, res2)

MAXITER = 100
betahat = beta_pilot
Xtmp = cbind(rep(1, N), X)
for (t in 1:MAXITER){
    resid = Y - Xtmp %*% betahat

    cur_gvals = evalScore(resid, knots, knot_gvals, knot_gprimes)
    cur_gprimes = evalScoreDeriv(resid, knots, knot_gvals, knot_gprimes)

    grad = -t(Xtmp) %*% cur_gvals
    hess = t(Xtmp) %*% diag(cur_gprimes) %*% Xtmp
    update = solve(hess, grad)

    betahat = betahat - 0.5 * update

    if (sum(grad^2) < 1e-7){
        print(sprintf("Final grad norm-sq %.5f", sum(grad^2)))
        break
    }
}