### Gaussian Process
setwd("~/Documents/SDS 383D/Exercise 3")
library(mvtnorm)
source('PointWiseCIPoly.R')

m0GP <- function(inputData,hyperB,hyperTau1sq,hyperTau2sq,covFunc = CovSqExp){
  #Input
  # @inputData: vector of input data values
  # @hyperB: hyper parameter b for covariance kernel, positive single scalar parameter, suggestion [.1,10]
  # @hyperTau1sq: hyper parameter tau1sq for covariance kernel, positive single scalar parameter, suggestion [.01,1]
  # @hyperTau2sq: hyper parameter tau2sq for covariance kernel, positive single scalar parameter, often kept close to 0
  #Output
  # Returns a vector of values representing the random process evaluated
  #   at each inputData point
  
  ## Sets mean vector and initializes Covariance Matrix
  meanFx = rep(0,length(inputData))
  covFx = matrix(NA,nrow=length(inputData),ncol=length(inputData))

  ## Fills out Covariance Matrix
  for(i in 1:length(inputData)){
    for(j in 1:length(inputData)){
      covFx[i,j] = covFunc(inputData[i],inputData[j],hyperB,hyperTau1sq,hyperTau2sq)
    }
  }
  ## Samples multivariate normal distribution and returns values
  return(rmvnorm(1,meanFx,covFx)[1,])
}

CovSqExp <- function(x1,x2,b, tau1sq, tau2sq){
  #Input
  # @x1, single scalar data value
  # @x2, single scalar data value
  # @b, positive single scalar parameter, suggestion [.1,10]
  # @tau1sq, positive single scalar parameter, suggestion [.01,1]
  # @tau2sq, positive single scalar parameter, often kept close to 0
  #Output
  # Returns Squared Exponential Covariance, a positive value
  #   between 0 and tau1sq + tau2sq
  
  ## Finds distance between data values
  dist = abs(x1-x2)
  
  ## Evaluates exponential term and kronecker delta function
  expTerm = exp(-1/2 * (dist/b)^2)
  delta = kroneckerDelta(x1,x2)
  
  ## Evaluates each term in sum separately
  firstTerm = tau1sq * expTerm
  secondTerm = tau2sq * delta
  
  ## Adds and returns
  return(firstTerm+secondTerm)
}

CovMat52 <- function(x1,x2,b,tau1sq,tau2sq){
  #Input
  # @x1, single scalar data value
  # @x2, single scalar data value
  # @b, positive single scalar parameter, suggestion [.1,10]
  # @tau1sq, positive single scalar parameter, suggestion [.01,1]
  # @tau2sq, positive single scalar parameter, often kept close to 0
  #Output
  # Returns Matern covariance of the two points using the parameter
  #   nu = 5/2, a positive value between 0 and tau1sq + tau2sq
  
  ## Finds distance between points
  dist = abs(x1-x2)
  
  ## Evaluates exponential term and taylor series term
  expTerm = exp(-sqrt(5)*dist/b)
  taylorTerm = 1 + sqrt(5)*dist/b + 5 * dist^2/(3*b^2)
  
  ## Evaluates kronecker delta function for data points
  delta = kroneckerDelta(x1,x2)
  
  ## valuates each term in sum separately
  term1 = tau1sq * taylorTerm * expTerm
  term2 = tau2sq * delta
  
  ## Sums and returns
  return(term1 + term2)
}

kroneckerDelta <- function(x1,x2){
  #Input: Two scalar values
  #Returns 1 if the input values are equal, otherwise returns 0
  if(x1 == x2) return(1)
  else return(0)
}

## Set scaling value, number of input points, and number of repititions to simulate process
c=1
reps=3
numPoints = 200

## Create Input Data
data = runif(numPoints)
data = data[order(data)]

## Set Hyperparameters
tau2sqd = 10^-6
tau1sqd = 1
b = .1

## Run Gaussian process for both covariance kernels
simData = replicate(reps,m0GP(data,b,tau1sqd,tau2sqd))
simData52 = replicate(reps,m0GP(data,b,tau1sqd,tau2sqd,CovMat52))

####################################################
#### Unneeded, but possibly useful later
## Find 95% confidence bands for both kernels
#bound95cov1 = sapply(1:numPoints,function(t) quantile(simData[t,],c(.025,.975)))
#bound95cov2 = sapply(1:numPoints,function(t) quantile(simData52[t,],c(.025,.975)))
####################################################

## Plot true values along with GP Approximation Bands
par(mfrow=c(1,2))
## Squared Exponential Covariance
plot(data,simData[,1],ylim=c(min(simData),max(simData)),xlim=c(0,1),main='Squared Exponential')
lines(data,simData[,1],col='red')
invisible(lapply(1:reps,function(t) points(data,simData[,t])))
invisible(lapply(1:reps,function(t) lines(data,simData[,t],col='red')))

## Matern Covariance
plot(data,simData52[,1],ylim=c(min(simData52),max(simData52)),xlim=c(0,1),main='Matern 5/2')
lines(data,simData52[,1],col='red')
invisible(lapply(1:reps,function(t) points(data,simData52[,t])))
invisible(lapply(1:reps,function(t) lines(data,simData52[,t],col='red')))


####################################################
## Altering Value of b
bseq = c(.01,.1,1,10)
par(mfrow=c(2,length(bseq)))
## Create Input Data
data = runif(numPoints)
data = data[order(data)]

for(b in bseq){
  ## Run Gaussian process for both covariance kernels
  simData = replicate(reps,m0GP(data,b,tau1sqd,tau2sqd))
  
  ## Plot true values along with GP Approximation Bands
  plot(data,simData[,1],ylim=c(min(simData),max(simData)),xlim=c(0,1),main='Squared Exponential',xlab=paste('b = ',toString(b)))
  lines(data,simData[,1],col='red')
  invisible(lapply(1:reps,function(t) points(data,simData[,t])))
  invisible(lapply(1:reps,function(t) lines(data,simData[,t],col='red')))
}

for(b in bseq){
  ## Run Gaussian process for both covariance kernels
  simData52 = replicate(reps,m0GP(data,b,tau1sqd,tau2sqd,CovMat52))
  
  ## Plot true values along with GP Approximation Bands
  plot(data,simData52[,1],ylim=c(min(simData52),max(simData52)),xlim=c(0,1),main='Matern 5/2',xlab=paste('b = ',toString(b)))
  lines(data,simData52[,1],col='red')
  invisible(lapply(1:reps,function(t) points(data,simData52[,t])))
  invisible(lapply(1:reps,function(t) lines(data,simData52[,t],col='red')))
}
####################################################

b=.1
####################################################
## Altering Value of tau1
tau1sqdSeq = c(.001,.01,.1,1)
par(mfrow=c(2,length(tau1sqdSeq)))
## Create Input Data
data = runif(numPoints)
data = data[order(data)]

for(tau1sqd in tau1sqdSeq){
  ## Run Gaussian process for both covariance kernels
  simData = replicate(reps,m0GP(data,b,tau1sqd,tau2sqd))
  
  ## Plot true values along with GP Approximation Bands
  plot(data,simData[,1],ylim=c(-3,3),xlim=c(0,1),main='Squared Exponential',xlab=paste('Tau1 = ',toString(tau1sqd)))
  lines(data,simData[,1],col='red')
  invisible(lapply(1:reps,function(t) points(data,simData[,t])))
  invisible(lapply(1:reps,function(t) lines(data,simData[,t],col='red')))
}

for(tau1sqd in tau1sqdSeq){
  ## Run Gaussian process for both covariance kernels
  simData52 = replicate(reps,m0GP(data,b,tau1sqd,tau2sqd,CovMat52))
  
  ## Find 95% confidence bands for both kernels
  bound95cov1 = sapply(1:numPoints,function(t) quantile(simData52[t,],c(.025,.975)))
  
  ## Plot true values along with GP Approximation Bands
  plot(data,simData52[,1],ylim=c(-3,3),xlim=c(0,1),main='Matern 5/2',xlab=paste('Tau1 = ',toString(tau1sqd)))
  lines(data,simData52[,1],col='red')
  invisible(lapply(1:reps,function(t) points(data,simData52[,t])))
  invisible(lapply(1:reps,function(t) lines(data,simData52[,t],col='red')))
}
####################################################