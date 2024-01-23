library(ggplot2)
library(gridExtra)
library(gnorm) # to generate generalized Gaussian

#rm(list = ls())  # clear all variables in workspace
#graphics.off()  # clear all plots

N <- 800 # sample size of each fold
d <- 2 # dimension of X
beta_0 <- 3*runif(d) + 1
mu <- 2
lambda0 <- 0.0001 # penalization coefficient

## Generate the data
X <- rnorm(N*d, 0, 1)
X <- array(X, c(N, d))
Y <- X%*%beta_0 + mu
# Additive noise
# Noise <- runif(N,-1,1)
# Noise <- rnorm(N,0,1)
Noise <- rnorm(N, 0, 1)
negNoise <- (Noise < 0)
# Asymmetric log-concave noise with mean = median = 0
Noise <- runif(N, -2*sqrt(2/pi),0)*negNoise + Noise*(1-negNoise)
Y <- Y + Noise

##the normalized penalization coefficient
lambda <- lambda0*2*N*sd(X)^3 #sd(X)^3: to make lambda scale invariant

###use OLS as pilot
Xtemp <- cbind(rep(1,N),X)
betapilot <- solve(t(Xtemp)%*%Xtemp,t(Xtemp)%*%Y)#betapilot[1] is the intercept
Residualsorted <- sort(Y-Xtemp%*%betapilot)
deltaResidual <- diff(Residualsorted)

u <- c(0,rep(1,N))
L <- t(matrix(replicate(N,c(2,deltaResidual)),ncol=N))
L[upper.tri(L)]  <-  0
A <- diag(x=1, nrow = N, ncol = N+1)
diag(A[-1,3:(N+1)]) <- 1
A[1,1] <-  -1
D <- diag(1/deltaResidual*lambda) #lambda= lambda0*2*N*sd(X)^3
B <- diag(x=1, nrow = N-1, ncol = N)
diag(B[,-1]) <-  -1
#D and B above are not the B and D in the tex file;
#we only use them to calculate the non-zero lower right block of the Hessian matrix

Hess <- t(B)%*%D%*%B #non-zero lower right block of the Hessian matrix
Hess <- cbind(rep(0,N),Hess)
Hess <- rbind(rep(0,N+1),Hess)
Hess <- 0.25*crossprod(L%*%A)+Hess #Hess=0.5*Hessian of the objective function, (Hess%*%G-u)=0.5*gradient
##initialization
#Gtemp = (-g(T_1),g'(T_1),g'(T_2),...,g'(T_N))
# g is -psi, monotone increasing
# sum_i g*(Ti)=0, g*(T1)<0, -g*(T1)>0, this constraint is redundant but makes the code cleaner 
Gtemp <- solve(Hess,u) #if no shape constraint, this is the g*
Gtemp[Gtemp<0] <- 0 #projection
##projected Newton
###setting max number of iterations
M <- 9
for (j in 1:M) {
  #cat("Newton loop loss",t(Gtemp)%*%Hess%*%Gtemp-2*t(u)%*%Gtemp,"i=",i,"\n")
  alpha <- 1 #reset step size
  Free <- (Gtemp>0)|(Hess%*%Gtemp-u<=0)
  ###backtracking
  Gtemp1 <- Gtemp
  gradient <- (Hess%*%Gtemp-u)[Free]
  Update <- solve(Hess[Free,Free],gradient)
  Gtemp1[Free] <- Gtemp[Free]-alpha*Update
  Gtemp1[Free][Gtemp1[Free]<0] <- 0#projection
  # backtracking
  while (t(Gtemp1)%*%Hess%*%Gtemp1-2*t(u)%*%Gtemp1>t(Gtemp)%*%Hess%*%Gtemp-2*t(u)%*%Gtemp+0.5*t((Gtemp1-Gtemp)[Free])%*%(2*gradient)) {
    alpha <- alpha*0.7
    Gtemp1[Free] <- Gtemp[Free]-alpha*Update
    Gtemp1[Free][Gtemp1[Free]<0] <- 0#projection
  }
  ####termination criterion
  if (sum((Gtemp-Gtemp1)^2)<N*1e-13){
    Gtemp <- Gtemp1
    break
  }
  Gtemp <- Gtemp1
}
#######G= (-g*(T1),g*'(T1),...,g*'(Tn))#######
G <- Gtemp

#######gvalues = (g*(T1),...,g*(Tn))#######
gvalues <- G[-1]
gvalues <- gvalues[-1]+gvalues[-N]
gvalues <- 0.5*gvalues*deltaResidual
for (j in 2:(N-1)) {
  gvalues[j] <- gvalues[j]+gvalues[j-1]
}
gvalues <- c(-G[1],-G[1]+gvalues)
#sum(gvalues) #check ==0
#############plot g*############
tailgrid <- (1:(N/100))*mean(deltaResidual)
tailY <- gvalues[N]+tailgrid*G[N+1]
Negtailgrid <- (-(N/100):-1)*mean(deltaResidual)
NegtailY <- gvalues[1]+Negtailgrid*G[2]
Xplot <- c(Negtailgrid+Residualsorted[1],Residualsorted,tailgrid+Residualsorted[N])
Yplot <- c(NegtailY,gvalues,tailY)

mydf <- data.frame(X = Xplot, Y = Yplot)
plot0.obj <- ggplot(data=mydf)+geom_point(aes(x=X, y=Y), color="red")+geom_line(aes(x=X, y=(sign(X)+1)*X/2), color="blue")
grid.arrange(plot0.obj)
####################################
###gbeta= - \nabla_{beta} sum_i loglikelihd(Y_i-<X_i,beta>), a vector
gbeta <- function(betahat,X,Y,cuts,gvalues,GG,deltaResidual){
  #deltaResidual <- diff(cuts)
  X <- cbind(rep(1,length(Y)),X)
  residual <- Y-X%*%betahat
  delta <- (residual)
  index <- findInterval(delta, cuts)
  index[index==0] <- 1#left linear region
  Values <- gvalues[index] #g(left end point)
  distoleft <- delta-cuts[index]
  Values <- Values+distoleft*GG[index] #linear part
  quadpts <- index!=length(cuts)&index!=1
  quadindex <- index[quadpts]
  Values[quadpts] <- Values[quadpts]+0.5*(GG[quadindex+1]-GG[quadindex])*distoleft[quadpts]^2/deltaResidual[quadindex]
  Values <- t(Values)%*%(X)
  return(t(Values))
}
###gprimebeta = Hessian{ - sum_i loglikelihd(Y_i-<X_i,beta>)}
gprimebeta <- function(betahat,X,Y,cuts,gvalues,GG,deltaResidual){
  #deltaResidual <- diff(cuts)
  X <- cbind(rep(1,length(Y)),X)
  delta <- (Y-X%*%betahat)
  index <- findInterval(delta, cuts)
  index[index==0] <- 1#left linear region
  Values <- GG[index] #g(left end point)
  quadpts <- index!=length(cuts)&index!=1
  quadindex <- index[quadpts]
  distoleft <- delta-cuts[index]
  Values[quadpts] <- Values[quadpts]+(GG[quadindex+1]-GG[quadindex])*distoleft[quadpts]/deltaResidual[quadindex]
  Values <- -t(X)%*%diag(Values)%*%X
  return(Values)
}

betahat <- betapilot #initialization
for (l in 1:8) {
  gvalue1 <- gbeta(betahat,X,Y,Residualsorted,gvalues,G[-1],deltaResidual)
  Update <- solve(gprimebeta(betahat,X,Y,Residualsorted,gvalues,G[-1],deltaResidual),gvalue1)
  betahat <- betahat-Update
}
#print(sum(gbeta(betahat,X,Y,Residualsorted,gvalues,G[-1],deltaResidual)^2))

abs(betahat-c(mu,beta_0))
abs(betapilot-c(mu,beta_0))
