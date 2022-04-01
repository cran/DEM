#' The DMOEM is an overrelaxation algorithm in distributed manner, which is used to solve the parameter estimation of multivariate Gaussian mixture model.
#'
#' @param y is a data matrix
#' @param M is the number of subsets
#' @param seed is the recommended way to specify seeds
#' @param alpha0 is the initial value of the mixing weight under the EM algorithm
#' @param mu0 is the initial value of the mean under the EM algorithm
#' @param sigma0 is the initial value of the covariance under the EM algorithm
#' @param MOEMalpha0 is the initial value of the mixing weight under the MOEM algorithm
#' @param MOEMmu0 is the initial value of the mean under the MOEM algorithm
#' @param MOEMsigma0 is the initial value of the covariance under the MOEM algorithm
#' @param omega is the overrelaxation factor
#' @param i is the number of iterations
#' @param epsilon is the threshold value
#'
#' @return DMOEMalpha,DMOEMmu,DMOEMsigma,DMOEMtime
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
#' MOEMalpha0= alpha1
#' MOEMmu0=mu1
#' MOEMsigma0=sigma1
#' omega=0.15
#' i=10
#' epsilon=0.005
#' DMOEM(y,M,seed,alpha0,mu0,sigma0,MOEMalpha0,MOEMmu0,MOEMsigma0,omega,i,epsilon)

DMOEM=function(y,M,seed,alpha0,mu0,sigma0,MOEMalpha0,MOEMmu0,MOEMsigma0,omega,i,epsilon){
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
MOEMalpha=MOEMalpha0 
MOEMmu=MOEMmu0   
MOEMsigma=MOEMsigma0
den=matrix(rep(0, K*nm), nrow = nm) 
prob=matrix(rep(0, K*nm), nrow = nm) 
weight=matrix(rep(0, K*nm), nrow = nm) 
r= c(rep(0,K))
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
oldMOEMalpha=MOEMalpha
oldMOEMmu=MOEMmu
oldMOEMsigma=MOEMsigma
for (k in 1:K){
varmat = matrix(0, ncol=ncol(y), nrow=ncol(y)) 
for (j in 1:nm){
varmat = varmat + prob[j,k] *(y1[j,]-mu[k,])%*%t((y1[j,]-mu[k,])) 
}
alpha[k]=prob1[k]/nm 
mu[k,] = (t(y1) %*% prob[,k]) / prob1[k] 
sigma[[k]] = varmat/prob1[k]  
}
for (k in 1:K){
r[k]=alpha[k]/oldMOEMalpha[k] 
}
omegat= omega*(min(r))/(1+omega-omega*(min(r))) 
for (k in 1:K){
MOEMalpha[k]=(1+omegat)*alpha[k]-omegat*oldMOEMalpha[k]
MOEMmu[k,]=(1+omega)*mu[k,]-omega*oldMOEMmu[k,]  
MOEMsigma[[k]]=(1+omega)*sigma[[k]]-omega*oldMOEMsigma[[k]] 
}
Sigma=c(rep(0,K))
for (k in 1:K){
Sigma[k]=max(abs(MOEMsigma[[k]]-oldMOEMsigma[[k]]))
}
    if(max(abs(MOEMalpha-oldMOEMalpha))<epsilon &
       max(abs(MOEMmu-oldMOEMmu))<epsilon & 
       max(Sigma)<epsilon)break
cat(
   "step",step,"\n",
   "MOEMalpha",MOEMalpha,"\n",
   "MOEMmu",MOEMmu,"\n")
}
alphaM=MOEMalpha+alphaM  
muM=MOEMmu+muM 
for (k in 1:K){
sigmaM[[k]]=MOEMsigma[[k]]+sigmaM[[k]]
}
}
)
time=time1/M
alphamao=alphaM/M 
mumao=muM/M
sigmamao=list()
for (k in 1:K){
sigmamao[[k]]=sigmaM[[k]]/M
}
return(list(DMOEMalpha=alphamao,DMOEMmu=mumao,DMOEMsigma=sigmamao,DMOEMtime=time))
}