## Cheese
library(mvtnorm)
library(MCMCpack)
setwd("~/Documents/SDS 383D/Exercise 4")
data = read.csv('~/Documents/SDS 383D/Data/cheese.csv')

## Define Data Variables
lP = log(data$price)
lQ = log(data$vol)
P = data$price
Q = data$vol
D = data$disp

## Because I use python and hate typing out length every time
len <- function(x) {length(x)}

## Setup size variables, and separate datapoints by store
N = len(P)
storeList = unique(data$store)
n = len(storeList)
lPi = sapply(1:n,function(t) lP[data$store == storeList[t]])
lQi = sapply(1:n,function(t) lQ[data$store == storeList[t]])
Di = sapply(1:n,function(t) D[data$store == storeList[t]])
mi = sapply(lPi,len)

## Create Input Matrices X (list of matrices where each matrix is for a different store)
X = sapply(1:n,function(t) matrix(c(rep(1,mi[t]),lPi[[t]],Di[[t]],lPi[[t]]*Di[[t]]),ncol=4))
XFull = matrix(c(rep(1,N),lP,D,D*lP),ncol=4)

## Setup hyperparameters for overall variance
a = 1
b = 1

## Initialize parameters
Mu = matrix(c(0,-2,0,0),nrow=4)
Sigma = diag(4)
SigmaInv = solve(Sigma)
Beta_i = rmvnorm(1,Mu,Sigma)[1,]
Lambda_Q = rgamma(1,a,b)

## Create Storage variables to keep every iterated parameter value
Store = 100
muStore = array(NA,dim=c(Store,4))
BetaStore = array(NA,dim=c(Store,4,n))

## Loop through Gibbs framework
for(i in 1:Store){
  ## Set up beta parameters and sample betas
  Beta_Var = sapply(1:n, function(t) solve(Lambda_Q * t(X[[t]]) %*% X[[t]] + SigmaInv))
  Beta_Mean = sapply(1:n, function(t) matrix(Beta_Var[,t],nrow=4) %*% (Lambda_Q * t(X[[t]]) %*% lQi[[t]] + SigmaInv %*% Mu))
  Beta_i = sapply(1:n,function(t) rmvnorm(1,Beta_Mean[,t],matrix(Beta_Var[,t],nrow=4)))
  
  ## Sample Beta Mean (Mu)
  Mu = t(rmvnorm(1,rowMeans(Beta_i),Sigma))
  
  ## Setup Beta Covariance Matrix (Sigma) and Sample
  BMinMu = sapply(1:n,function(t) Beta_i[,t] - Mu)
  SigBMat = matrix(rowMeans(sapply(1:n,function(t) BMinMu[,t] %*% t(BMinMu[,t]))),nrow=4)
  Sigma = riwish(n/2 - 1, SigBMat)
  
  ## Setup Lambda Parameters
  lQMinXB = sapply(1:n,function(t) (lQi[[t]] - X[[t]] %*% matrix(Beta_i[,t],nrow=4))^2)
  lQMinXBSum = sapply(lQMinXB,sum)
  Lambda_Q = rgamma(1,a + N/2,b + 1/2*sum(lQMinXBSum)/N)
  
  ## Store sampled values
  muStore[i,] = Mu
  BetaStore[i,,] = Beta_i
}

pdf(file='Plots')  

## Plot Posterior Estimates of Betas for store 1
par(mfrow=c(2,2))
hist(BetaStore[,1,1],breaks = 20)
hist(BetaStore[,2,1],breaks = 20)
hist(BetaStore[,3,1],breaks = 20)
hist(BetaStore[,4,1],breaks = 20)

## Plot Posterior Average (over all stores) of Betas 
par(mfrow=c(2,2))
postBetaIs = colMeans(BetaStore) 
hist(postBetaIs[1,],breaks = 20)
hist(postBetaIs[2,],breaks = 20)
hist(postBetaIs[3,],breaks = 20)
hist(postBetaIs[4,],breaks = 20)

## Plot Sampled Mean Parameters used for Betas
par(mfrow=c(2,2))
hist(muStore[,1],breaks=20)
hist(muStore[,2],breaks=20)
hist(muStore[,3],breaks=20)
hist(muStore[,4],breaks=20)

## Calculate Posterior Estimates of Betas agnostic of store
postBetas = rowMeans(postBetaIs)

## Show relationship between Betas by plotting Mus (Shows correlations)
par(mfrow=c(2,3))
for(i in 1:3){
  for(j in (i+1):4){
    plot(muStore[,i],muStore[,j],main=paste('Scatterplot of Mu',toString(i),'and Mu',toString(j)))
  }
}

## Creates estimates for every point using posterior parameters
Yhati = sapply(1:n, function(t) X[[t]] %*% postBetaIs[,t])

## Simple helper function to make display function determine color in plots
colorTrans <- function(vec01){
  colors = rep(NA,len(vec01))
  i = 1
  for(bit in vec01){
    if(bit) colors[i] = 'red'
    else colors[i] = 'blue'
    i = i + 1
  }
  return(colors)
}

## Get colored up
colors = sapply(Di,colorTrans)

## Plot store by store regressions alongside true data
par(mfrow=c(1,1))
for(i in 1:n) {
  ## Plot True Data
  plot(lPi[[i]],lQi[[i]],col=colors[[i]],main=storeList[i],xlab='Log Price',ylab='Log Volume')
  
  ## Separate data by display
  ixOn = Di[[i]]
  pOn = lPi[[i]][ixOn == T]
  pOff = lPi[[i]][ixOn == F]
  yOn = Yhati[[i]][ixOn == T]
  yOff = Yhati[[i]][ixOn == F]
  
  ## Plot regression curves
  curve(postBetaIs[1,i] + postBetaIs[3,i] + x * (postBetaIs[2,i] + postBetaIs[4,i]),add=T,col='red')
  curve(postBetaIs[1,i] + x * (postBetaIs[2,i]),add=T,col='blue')
}

## Calculate optimal price values
c=1
par(mfrow=c(1,2))
newPNonDisp = c * postBetaIs[2,] / (postBetaIs[2,] + 1)
newPDisp = c* (postBetaIs[2,] + postBetaIs[4,])/(postBetaIs[2,] + postBetaIs[4,] + 1)
## Some values can be weird, but in general the evalues are reasonable
hist(newPNonDisp,breaks=40,main='Histogram of Optimal Non-Display Prices')
hist(newPDisp,breaks=40,main='Histogram of Optimal Display Prices')
median(newPNonDisp)
median(newPDisp)
dev.off()

