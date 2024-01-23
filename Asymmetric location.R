#library(RSpectra) #svds() used in gradient descent
library(ggplot2)
library(gridExtra)
library(gnorm) #to generate generalized Gaussian
#library(quadprog) #solve.QP(), off-the-shelf quadratic programming

#rm(list = ls())  # clear all variables in workspace
#graphics.off()  # clear all plots

##set the penalization coefficient
lambda0 <-0.0001
##sample size
N<-1000

##generate data
# X<-rexp(N)
# X<-X*sign(runif(N,-1,1))
# X<-runif(N,-1,1)
X<-rnorm(N,0,1)
# X<-rgnorm(N, mu = 0, alpha = 1, beta = 9) #generalized Gaussian
# X<-runif(N/2,-5,-3)
# X<-c(X,runif(N/2,2.5,5.5))
#log-concave, mean=median=0, asymmetric
# X<-rnorm(N,0,1)
# negX<-X<0
# X<-runif(N,-2*sqrt(2/pi),0)*negX+X*(1-negX)

lambda<-lambda0*2*N*sd(X)^3 #sd(X)^3: to make lambda scale invariant
##setting thetapilot
thetapilot<-mean(X)
#thetapilot<-0 #oracle if the noise is symmetric

Xsorted<-sort(X-thetapilot)#T_1,T_2,...,T_N
deltaX<-diff(Xsorted)#deltaX[i]=T_{i+1}-T_i


#G=(-g(T_1),g'(T_1),g'(T_2),...,g'(T_N))
# g is -psi, monotone increasing
##small scale projected Newton
u<-c(0,rep(1,N))
L<-t(matrix(replicate(N,c(2,deltaX)),ncol=N))
L[upper.tri(L)] <- 0
A<-diag(x=1, nrow = N, ncol = N+1)
diag(A[-1,3:(N+1)])<-1
A[1,1]<--1 #reparametrization, G[1]=-g(T_1)>0

####test if the matrix expression is correct####
# G<-runif(N+1, min=-0.3, max=0.2)
# GG<-G[-1]
# GG<-GG[-1]+GG[-N]
# GG<-0.5*GG*deltaX
# for (j in 2:(N-1)) {
#   GG[j]<-GG[j]+GG[j-1]
# }
# GG<-c(-G[1],GG-G[1])
# sum(GG^2)-t(G)%*%crossprod(L%*%A)%*%G*0.25
################################################

D<-diag(1/deltaX*lambda) #lambda=lambda0*2*N*sd(X)^3
B<-diag(x=1, nrow = N-1, ncol = N)
diag(B[,-1])<--1
#D and B above are not the B and D in the tex file;
#we only use them to calculate the non-zero lower right block of the Hessian matrix

####test if the matrix expression is correct####
# G<-runif(N+1, min=-0.3, max=0.2)
# GG<-G[-1]
# GG<-GG[-1]-GG[-N]
# GG<-GG^2/deltaX
# Hess<-t(B)%*%D%*%B #Hess=0.5Hessian, (Hess%*%G-u)=0.5gradient
# Hess<-cbind(rep(0,N),Hess)
# Hess<-rbind(rep(0,N+1),Hess)
# lambda*sum(GG)-t(G)%*%Hess%*%G
################################################

Hess<-t(B)%*%D%*%B #non-zero lower right block of the Hessian matrix
Hess<-cbind(rep(0,N),Hess)
Hess<-rbind(rep(0,N+1),Hess)
Hess<-0.25*crossprod(L%*%A)+Hess #Hess=0.5*Hessian of the objective function, (Hess%*%G-u)=0.5*gradient
#Scaling<-solve(Hess)

##initialization, G=(-g(T_1),g'(T_1),g'(T_2),...,g'(T_N))
# g is -psi, monotone increasing
G<-solve(Hess,u)
G[G<0]<-0 #projection
###projected Newton
M<-20
for (i in 1:M) {
  cat("Newton loop loss",t(G)%*%Hess%*%G-2*t(u)%*%G,"\n")
  alpha<-1 #reset step size; more aggressive
  # Fixed<-(G==0)*(Hess%*%G-u>0)
  # Free<-!Fixed
  Free<-(G>0)|(Hess%*%G-u<=0)
  
  ###no backtracking for alpha
  # G[Free]<-G[Free]-alpha*Scaling[Free,Free]%*%(Hess%*%G-u)[Free]
  # G[Free][G[Free]<0]<-0#projection
  
  ###backtracking
  Gtemp<-G
  gradient<-(Hess%*%G-u)[Free]
  Update<-solve(Hess[Free,Free],gradient)
  Gtemp[Free]<-G[Free]-alpha*Update
  Gtemp[Free][Gtemp[Free]<0]<-0#projection
  while (t(Gtemp)%*%Hess%*%Gtemp-2*t(u)%*%Gtemp>t(G)%*%Hess%*%G-2*t(u)%*%G+0.5*t((Gtemp-G)[Free])%*%(2*gradient)) {
    alpha<-alpha*0.8
    Gtemp[Free]<-G[Free]-alpha*Update
    Gtemp[Free][Gtemp[Free]<0]<-0#projection
  }
  if(sum((Gtemp-G)^2)<1e-30){
    cat("the last update has Euclidean norm",sum((Gtemp-G)^2))
    G<-Gtemp
    break
    }
  G<-Gtemp
}


###########off-the-shelf QP######################
#same result
#qp <- solve.QP(Dmat=Hess, dvec=u, Amat=diag(x=1, nrow = N+1, ncol = N+1), bvec=rep(0,N+1), meq = 0)
#sum((qp$solution-G)^2)
############################################

######## g* ########
gvalues<-G[-1]
gvalues<-gvalues[-1]+gvalues[-N]
gvalues<-0.5*gvalues*deltaX
for (j in 2:(N-1)) {
  gvalues[j]<-gvalues[j]+gvalues[j-1]
}
gvalues<-c(-G[1],-G[1]+gvalues) #gvalues = (g*(T1),...,g*(Tn))
########d/dx G = g*, Negloglikelihd = (G(T1),...,G(Tn))
GG<-G[-1]
Negloglikelihd<-gvalues[-N]*deltaX+(GG[-N]/2+(GG[-1]-GG[-N])/6)*deltaX^2
for (j in 2:(N-1)) {
  Negloglikelihd[j]<-Negloglikelihd[j]+Negloglikelihd[j-1]
}
normalization<-0 #normalization[i] = log(integral exp(-G)) where G(T_{N/2})=0
Negloglikelihd<-c(0,Negloglikelihd)
Negloglikelihd<-Negloglikelihd-Negloglikelihd[N/2]+normalization

#############PLOT############
tailgrid<-(1:(N/100))*mean(deltaX)
Negtailgrid<-(-(N/100):-1)*mean(deltaX)
Xplot<-c(Negtailgrid+Xsorted[1],Xsorted,tailgrid+Xsorted[N])
############# g ############
##the linear tail
NStailY<-gvalues[N]+tailgrid*G[N+1]
NSNegtailY<-gvalues[1]+Negtailgrid*G[2]
NSYplot<-c(NSNegtailY,gvalues,NStailY)
#############negative log-likelihood############
NLLtailY<-Negloglikelihd[N]+tailgrid*gvalues[N]+tailgrid^2/2*G[N+1]
NLLNegtailY<-Negloglikelihd[1]+Negtailgrid*gvalues[1]+Negtailgrid^2/2*G[2]
NLLYplot<-c(NLLNegtailY,Negloglikelihd,NLLtailY)
Xmin<-min(Xsorted)
Xmax<-max(Xsorted)

#mydf = data.frame(X = Xplot[(N/4):(N*3/4)],NScore = NSYplot[(N/4):(N*3/4)],NLL = NLLYplot[(N/4):(N*3/4)])
mydf = data.frame(X = Xplot,NScore = NSYplot,NLL = NLLYplot)
plotNS.obj <- ggplot(data=mydf)+geom_line(aes(x=X, y=NScore), color="red")+
  geom_line(aes(x=X, y=X), color="blue")+ylab("Negative score function")+
  geom_vline(xintercept = Xmin,linetype="dotted")+ geom_vline(xintercept = Xmax,linetype="dotted")
plotNLL.obj <- ggplot(data=mydf)+geom_line(aes(x=X, y=NLL), color="red")+geom_line(aes(x=X, y=X^2/2), color="blue")+
  ylab("Negative log-likelihood")+ geom_vline(xintercept = Xmin,linetype="dotted")+ geom_vline(xintercept = Xmax,linetype="dotted")

# plotNS.obj <- ggplot(data=mydf)+geom_line(aes(x=X, y=NScore), color="red")+
#   geom_line(aes(x=X, y=(sign(X-thetapilot)+1)*(X-thetapilot)/2), color="blue")+ylab("Negative score function")+
#   geom_vline(xintercept = Xmin,linetype="dotted")+ geom_vline(xintercept = Xmax,linetype="dotted")
# plotNLL.obj <- ggplot(data=mydf)+geom_line(aes(x=X, y=NLL), color="red")+geom_line(aes(x=X, y=(sign(X-thetapilot)+1)*(X-thetapilot)^2/4), color="blue")+
#   ylab("Negative log-likelihood")+ geom_vline(xintercept = Xmin,linetype="dotted")+ geom_vline(xintercept = Xmax,linetype="dotted")

grid.arrange(plotNS.obj, plotNLL.obj, ncol=2)
####################################

########knots=(T1,...,Tn),gvalues=(g(T1),...,g(Tn)),GG=(g'(T1),...,g'(Tn))=G[-1],deltaX=(T2-T1,...,Tn-T{n-1})########
gtheta<-function(theta,X1,knots,gvalues,GG,deltaX){
  #deltaX<-diff(knots)
  delta<-(X1-theta)
  index<-findInterval(delta, knots)
  index[index==0]<-1 #left linear region
  distoleft<-delta-knots[index]
  Values<-gvalues[index] #g(left end point)
  Values<-Values+distoleft*GG[index] #linear part
  quadpts<-index!=length(knots)&index!=1
  quadindex<-index[quadpts]
  Values[quadpts]<-Values[quadpts]+0.5*(GG[quadindex+1]-GG[quadindex])*distoleft[quadpts]^2/deltaX[quadindex]
  return(sum(Values))
}
###gprimetheta is d/dx g=-d/dtheta g
gprimetheta<-function(theta,X1,knots,gvalues,GG,deltaX){
  #deltaX<-diff(knots)
  delta<-(X1-theta)
  index<-findInterval(delta, knots)
  index[index==0]<-1#left linear region
  Values<-GG[index] #g(left end point)
  quadpts<-index!=length(knots)&index!=1
  quadindex<-index[quadpts]
  distoleft<-delta-knots[index]
  Values[quadpts]<-Values[quadpts]+(GG[quadindex+1]-GG[quadindex])*distoleft[quadpts]/deltaX[quadindex]
  return(sum(Values))
}
# #######compare with the sample mean
# Thetahats<-NULL
# Thetatildes<-NULL
# for (j in 1:100) {
#   #X1<-runif(N,-1,1)
#   #X1<-rnorm(N)
#   #X1<-rgnorm(N, mu = 0, alpha = 1, beta = 4)
# 
#   X1<-abs(rnorm(N/2,0,1))
#   X1<-c(X1,runif(N/2,-2*sqrt(2/pi),0))
# 
#   thetatilde<-mean(X1)
#   thetahat<-thetatilde #initialization
#   while (abs(gtheta(thetahat,X1,Xsorted,gvalues,G[-1],deltaX))>1e-13*N) {
#     # #line search for step size
#     # alpha<-1
#     # Update<-gtheta(thetahat,X1,Xsorted,gvalues,G[-1],deltaX)/gprimetheta(thetahat,X1,Xsorted,gvalues,G[-1],deltaX)
#     # thetatemp<-thetahat+alpha*Update
#     # while (abs(gtheta(thetatemp,X1,Xsorted,gvalues,G[-1],deltaX))>abs(gtheta(thetahat,X1,Xsorted,gvalues,G[-1],deltaX))) {
#     #   alpha<-0.8*alpha
#     #   thetatemp<-thetahat+alpha*Update
#     # }
#     # thetahat<-thetatemp
#     ###no backtracking
#     Update<-gtheta(thetahat,X1,Xsorted,gvalues,G[-1],deltaX)/gprimetheta(thetahat,X1,Xsorted,gvalues,G[-1],deltaX)
#     thetahat<-thetahat+alpha*Update
#     if(abs(Update)<1e-12){break}
#   }
#   #print(abs(gtheta(thetahat,X1,Xsorted,gvalues,G[-1],deltaX)))
#   Thetahats<-c(Thetahats,abs(thetahat))
#   Thetatildes<-c(Thetatildes,abs(thetatilde)/sqrt(2))
#   #divide the |mean-theta0| by sqrt(2) to balance the extra 1000 samples used for score matching
# }
# mean(Thetahats)
# mean(Thetatildes)
