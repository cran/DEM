#' The DEM2 algorithm is a one-step average algorithm in distributed manner, which is used to solve the parameter estimation of multivariate Gaussian mixture model.
#'
#' @param y is a data matrix
#' @param M is the number of subsets
#' @param seed is the recommended way to specify seeds
#' @param alpha0 is the initial value of the mixing weight
#' @param mu0 is the initial value of the mean
#' @param sigma0 is the initial value of the covariance
#' @param i is the number of iterations
#' @param epsilon is the threshold value
#'
#' @return DEM2alpha,DEM2mu,DEM2sigma,DEM2time
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
#' M=5
#' seed=123
#' alpha0= alpha1
#' mu0=mu1
#' sigma0=sigma1
#' i=10
#' epsilon=0.005
#' DEM2(y,M,seed,alpha0,mu0,sigma0,i,epsilon)

DEM2=function(y,M,seed,alpha0,mu0,sigma0,i,epsilon){
n=nrow(y)
p=ncol(y)
K=length(alpha0)
nm=n/M 
alphaM=c(rep(0,K)) 
muM=matrix(rep(0, K*p), nrow = K) 
sigmaM=list()
for (k in 1:K){
sigmaM[[k]]=matrix(rep(0, p*p), nrow = p)
}
set.seed(seed)
mr=matrix(sample(c(1:n),n,replace=FALSE),nrow = M,ncol=nm,byrow=TRUE)
time1=system.time(for (m in 1:M) {
y1=y[mr[m,],] 
alpha=alpha0 
mu=mu0 
sigma=sigma0
den=matrix(rep(0, K*nm), nrow = nm) 
prob=matrix(rep(0, K*nm), nrow = nm) 
weight=matrix(rep(0, K*nm), nrow = nm) 
for (step in 1:i){
for (k in 1:K){
den[, k]=dmvnorm(y1, mu[k,], sigma[[k]], log=FALSE) 
weight[, k]=alpha[k] * den[, k] 
}
prob=weight/rowSums(weight)
prob1=colSums(prob) 
oldalpha=alpha
oldmu=mu
oldsigma=sigma
for (k in 1:K){
varmat = matrix(0, ncol=ncol(y), nrow=ncol(y)) 
for (j in 1:nm){
varmat = varmat + prob[j,k] *(y1[j,]-mu[k,])%*%t((y1[j,]-mu[k,])) 
}
alpha[k]=prob1[k]/nm 
mu[k,] = (t(y1) %*% prob[,k]) / prob1[k] 
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
alphaM=alpha+alphaM  
muM=mu+muM 
for (k in 1:K){
sigmaM[[k]]=sigma[[k]]+sigmaM[[k]]
}
}
)
alphamao=alphaM/M 
mumao=muM/M 
sigmamao=list()
for (k in 1:K){
sigmamao[[k]]=sigmaM[[k]]/M
}
time=time1/M
return(list(DEM2alpha=alphamao,DEM2mu=mumao, DEM2sigma=sigmamao,DEM2time=time))
}