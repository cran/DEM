#' The EM algorithm is used to solve the parameter estimation of multivariate Gaussian mixture model.
#'
#' @param y is a data matrix
#' @param alpha0 is the initial value of the mixing weight
#' @param mu0 is the initial value of the mean
#' @param sigma0 is the initial value of the covariance
#' @param i is the number of iterations
#' @param epsilon is the threshold value
#'
#' @return EMalpha,EMmu,EMsigma,EMtime
#' @export
#'

#' @examples
#' library(mvtnorm)
#' alpha1= c(rep(1/4,4)) 
#' mu1=matrix(0,nrow=4,ncol=4) 
#' for (k in 1:4){
#' mu1[4,]=c(runif(4,(k-1)*3,k*3)) 
#' }
#' sigma1=list()
#' for (k in 1:4){
#' sigma1[[k]]= diag(4)*0.1
#' }
#' y= matrix(0,nrow=200,ncol=4) 
#' for(k in 1:4){
#' y[c(((k-1)*200/4+1):(k*200/4)),] = rmvnorm(200/4,mu1[k,],sigma1[[k]]) 
#' }
#' alpha0= alpha1
#' mu0=mu1
#' sigma0=sigma1
#' i=10
#' epsilon=0.005
#' EM(y,alpha0,mu0,sigma0,i,epsilon)

EM=function(y,alpha0,mu0,sigma0,i,epsilon){
n=nrow(y)
p=ncol(y)
K=length(alpha0)
alpha=alpha0 
mu=mu0 
sigma=sigma0
time=system.time(for (step in 1:i){
den=matrix(rep(0, K*n), nrow = n) 
prob=matrix(rep(0, K*n), nrow = n) 
weight=matrix(rep(0, K*n), nrow = n) 
for (k in 1:K){
den[, k]=dmvnorm(y, mu[k,], sigma[[k]], log=F) 
weight[, k]=alpha[k] * den[, k] 
}
prob=weight/rowSums(weight)
prob1=colSums(prob)
oldalpha=alpha
oldmu=mu
oldsigma=sigma
for (k in 1:K){
varmat = matrix(0, ncol=ncol(y), nrow=ncol(y)) 
for (j in 1:n){
varmat = varmat + prob[j,k] *(y[j,]-mu[k,])%*%t((y[j,]-mu[k,])) 
}
alpha[k]=prob1[k]/n 
mu[k,] = (t(y) %*% prob[,k]) / prob1[k] 
sigma[[k]] = varmat/prob1[k]  
}
Sigma=c(rep(0,K))
for (k in 1:K){
Sigma[k]=max(abs(sigma[[k]]-oldsigma[[k]]))
}
if(max(abs(alpha-oldalpha))<epsilon &
       max(abs(mu-oldmu))<epsilon & 
       max(Sigma)<epsilon)break 
  cat(
   "step",step,"\n",
   "alpha",alpha,"\n",
   "mu",mu,"\n"
)
}
)
return(list(EMalpha=alpha, EMmu=mu, EMsigma=sigma,EMtime=time))
}
