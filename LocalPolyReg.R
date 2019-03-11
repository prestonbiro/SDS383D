#Local Polynomial Regression
library(emdbook)
setwd("~/Documents/SDS 383D")
data = read.csv('Data/utilities.csv')

gaussianKernel <- function(distance,bandwidth){
  #Input:
  #distance is a real number (can be negative, but will take magnitude)
  #bandwidth is a positive real number 
  #Output:
  #The distance divided by the bandwidth evaluated by the Gaussian kernel function
  #   K(x) = 1/sqrt(2 pi) * exp(-x^2/2)
  #eval is the nonnegative real value returned for the evaluated value
  dist = abs(distance)
  xEval = dist/bandwidth
  eval = 1/sqrt(2*pi) * exp(-xEval^2/2)
  return(eval)
}

localLinearSingle <- function(dataX,dataY,center,kernelFunc,bandH){
  #Input:
  # @data:
  # @center:
  # @coefVec:
  
  #Output:
  # @

  N = length(dataX)
  
  Rx = matrix(1,nrow=N,ncol=2)
  Rx[,2] = dataX - center
  W = diag(kernelFunc(dataX - center,bandH))
  ahat = solve(t(Rx) %*% W %*% Rx) %*% t(Rx) %*% W %*% dataY
  return(ahat[1])
}

localLinear <- function(totalData,dataY,kernelFunc,bandH){
  #Inputs/Output
  indeces = order(totalData)
  
  estVec = sapply(1:length(totalData),function(t) localLinearSingle(totalData[-t],dataY[-t],totalData[t],kernelFunc,bandH))
  plot(totalData,dataY)
  points(totalData,estVec,cex=.5,col='blue')
  lines(totalData[indeces],estVec[indeces],col='red')
  
  return(dataY - estVec)
}

localPolyDSingle <- function(dataX,dataY,center,kernelFunc,bandH,D){
  #Input:
  # @data:
  # @center:
  # @coefVec:
  
  #Output:
  # @
  
  N = length(dataX)
  Rx = matrix(1,nrow=N,ncol=D+1)
  for(i in 1:(D)){
    Rx[,i+1] = (dataX - center)^i
  }
  print(Rx)
  W = diag(kernelFunc(dataX - center,bandH))
  ahat = solve(t(Rx) %*% W %*% Rx) %*% t(Rx) %*% W %*% dataY
  return(ahat[1])
}

localPoly <- function(totalData,dataY,kernelFunc,bandH,D){
  #Inputs/Output
  indeces = order(totalData)
  
  estVec = sapply(1:length(totalData),function(t) localPolyDSingle(totalData[-t],dataY[-t],totalData[t],kernelFunc,bandH,D))
  print(estVec)
  plot(totalData,dataY)
  points(totalData,estVec,cex=.5,col='blue')
  lines(totalData[indeces],estVec[indeces],col='red')
  return(estVec)
}

par(mfrow=c(1,1))
X = data$temp
Y = data$gasbill/data$billingdays
a = localLinear(X,Y,gaussianKernel,2)
# b = localPoly(X,Y,gaussianKernel,2,1)
# c = localPoly(X,Y,gaussianKernel,2,2)
# d = localPoly(X,Y,gaussianKernel,2,3)

## Heteroscadicity Check
plot(X,a,main='Residual Plot',ylab='Residual Value')


