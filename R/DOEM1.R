#' The DOEM1 algorithm is an online EM algorithm in distributed manner, which is used to solve the parameter estimation of multivariate Gaussian mixture model.
#'
#' @param y is a data matrix
#' @param M is the number of subsets
#' @param seed is the recommended way to specify seeds
#' @param alpha0 is the initial value of the mixing weight
#' @param mu0 is the initial value of the mean
#' @param sigma0 is the initial value of the covariance
#' @param i is the number of iterations
#' @param epsilon is the threshold value
#' @param a represents the power of the reciprocal of the step size
#' @param b indicates that the M-step is not implemented for the first b data points
#' @param c represents online iteration starting at 1/c of the total sample size
#'
#' @return DOEM1alpha,DOEM1mu,DOEM1sigma,DOEM1time
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
#' M=2
#' seed=123
#' alpha0= alpha1
#' mu0=mu1
#' sigma0=sigma1
#' i=10
#' epsilon=0.005
#' a=1
#' b=10
#' c=2
#' DOEM1(y,M,seed,alpha0,mu0,sigma0,i,epsilon,a,b,c)

DOEM1=function(y,M,seed,alpha0,mu0,sigma0,i,epsilon,a,b,c){
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
mr=matrix(sample(1:n,n,replace=FALSE),nrow = M,ncol=nm,byrow=TRUE)
time1=system.time(for (m in 1:M) {
y1=y[mr[m,],]
n1=nm/c
y2=y1[c(1:n1),]
alpha=alpha0 
mu=mu0 
sigma=sigma0
den=matrix(rep(0, K*n1), nrow = n1) 
weight=matrix(rep(0, K*n1), nrow = n1) 
for (step in 1:i){
for (k in 1:K){
den[, k]=dmvnorm(y2, mu[k,], sigma[[k]], log=FALSE) 
weight[, k]=alpha[k] * den[, k] 
}
prob=weight/rowSums(weight)
prob1=colSums(prob) 
oldalpha=alpha
oldmu=mu
oldsigma=sigma
for (k in 1:K){
varmat = matrix(0, ncol=ncol(y), nrow=ncol(y)) 
for (j in 1:n1){
varmat = varmat + prob[j,k] *(y2[j,]-mu[k,])%*%t((y2[j,]-mu[k,])) 
}
alpha[k]=prob1[k]/n1
mu[k,] = (t(y2) %*% prob[,k]) / prob1[k] 
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
y3=y1[c((n1+1):nm),]
n2=nm-n1
den1=c(rep(0,K)) 
weight1=c(rep(0,K))
for (k in 1:K){
den1[k]=dmvnorm(t(y3[1,]), mu[k,], sigma[[k]], log=FALSE)
weight1[k]=alpha[k] * den1[k]
}
prob2=weight1/sum(weight1)
S1=c(rep(0,K)) 
S2=matrix(0,nrow=K,ncol=p) 
S3=list()
for (k in 1:K){
S3[[k]]=matrix(0,p,p)
}
for (k in 1:K){
S1[k]= prob2[k]
S2[k,]=prob2[k]*y3[1,]
S3[[k]]=prob2[k]*(y3[1,]%*%t(y3[1,]))
}
den2=c(rep(0,K)) 
weight2=c(rep(0,K))

for (t in 0:(n2-1)) {
for (k in 1:K){
den2[k]=dmvnorm(t(y3[t+1,]), mu[k,], sigma[[k]], log=FALSE)
weight2[k]=alpha[k] * den2[k]
}
prob3=weight2/sum(weight2)
gamma=1/(t+1)^a 
oldS1=S1
oldS2=S2
oldS3=S3
oldalpha=alpha
oldmu=mu
oldsigma=sigma
for (k in 1:K){
S1[k]=(1-gamma)*oldS1[k]+gamma*prob3[k]
S2[k,]=(1-gamma)*oldS2[k,]+gamma*(prob3[k]*y3[t+1,])
S3[[k]]=(1-gamma)*oldS3[[k]]+gamma*(prob3[k]*(y3[t+1,]%*%t(y3[t+1,])))
if(t<b){alpha[k]=alpha[k]
mu[k,]=mu[k,] 
sigma[[k]]=sigma[[k]]} 
else {
alpha[k]=S1[k] 
mu[k,]=S2[k,]/S1[k] 
sigma[[k]]=S3[[k]]/S1[k]-S2[k,]%*%t(S2[k,])/(S1[k]^2)} 
}
cat(
   "alpha",alpha,"\n",
   "mu",mu,"\n")
}
alphaM=alpha+alphaM  
muM=mu+muM 
for (k in 1:K){
sigmaM[[k]]=sigma[[k]]+sigmaM[[k]]
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
return(list(DOEM1alpha=alphamao, DOEM1mu=mumao, DOEM1sigma=sigmamao,DOEM1time=time))
}