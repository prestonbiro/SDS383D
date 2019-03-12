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

localLinearSingle <- function(dataX,dataY,center,bandH,kernelFunc = gaussianKernel){
  #Input:
  # @dataX: a vector of real input data values
  # @dataY: a vector of real output data values
  # @center: single real value which to evaluate the smoothing
  # @bandH: bandwidth value, positive
  # @kernelFunc: kernel function used to evaluate the smoothing weights, defaulted to standard normal
  #Output:
  # @Returns the predicted y value at the 'center' point given after smoothing is conducted, single real value
  
  w = wix(dataX,center,bandH,kernelFunc)
  return(sum(w * dataY))
}

wix <- function(data,center,h,kernelF = gaussianKernel){
  #Input:
  # @data: a vector of real input data values
  # @center: single real value which to evaluate the smoothing
  # @h: bandwidth value, positive
  # @kernelFunc: kernel function used to evaluate the smoothing weights, defaulted to standard normal
  #Output:
  # @Returns the normalized weight vector for an individual x value (center) given the other x values (data)

  s1 = sjx(data,center,h,1)
  s2 = sjx(data,center,h,2)
  kernelEvals = sapply((center-data),function(t) kernelF(t,h))
  w = kernelEvals*(s2 - (data - center)*s1)
  return(w/sum(w))
}

sjx <- function(data,center,h,j,kernelF = gaussianKernel){
  #Input:
  # @data: a vector of real input data values
  # @center: single real value which to evaluate the smoothing
  # @h: bandwidth value, positive
  # @j: single integer valued exponent of (data - center) part of equation
  # @kernelFunc: kernel function used to evaluate the smoothing weights, defaulted to standard normal
  #Output:
  # @Helper function for wix, evaluating to a real value (possibly always positive?)
  
  kernelEvals = sapply((center-data),function(t) kernelF(t,h))
  distJ = (data - center)^j
  return(sum(kernelEvals*distJ))
}

localLinear <- function(dataX,dataY,bandH,kernelFunc = gaussianKernel){
  #Input:
  # @dataX: a vector of real input data values
  # @dataY: a vector of real output data values
  # @bandH: bandwidth value, positive
  # @kernelFunc: kernel function used to evaluate the smoothing weights, defaulted to standard normal
  #Output:
  # @Returns the predicted y values at each point using local smoothing is conducted, vector of real values
  
  ## Find predicted values
  estVec = sapply(1:length(dataX),function(t) localLinearSingle(dataX,dataY,dataX[t],bandH,kernelFunc))
  
  ## Plot true data with predicted values overlaid
  plot(dataX,dataY)
  points(dataX,estVec,cex=.5,col='blue')
  lines(dataX[indeces],estVec[indeces],col='red')
  
  return(estVec)
}

## Reset graphing and initialize data, also get indices for order X (helpful later)
par(mfrow=c(1,1))
X = data$temp
Y = data$gasbill/data$billingdays
indices = order(X)

## Get predicted values using local linear smoother
yhat = localLinear(X,Y,3)

## Optimizing Bandwidth
minH <- function(dataX,dataY,ker = gaussianKernel){
  #Input:
  # @dataX: a vector of real input data values
  # @dataY: a vector of real output data values
  # @kernelFunc: kernel function used to evaluate the smoothing weights, defaulted to standard normal
  #Output:
  # @Returns the value of bandwidth that minimizes the LOOCV value on a bandwidth range from .5 to 10.
  
  ## Initialize bandwidth options worth checking. Should alter based on data set
  hseq = seq(.5,10,.25)
  
  ## Find the LOOCV values for each value in hseq, find index of the minimum value
  sumRes = sapply(hseq,function(h) loocv(dataX,dataY,h,ker))
  ix = match(min(sumRes),sumRes)
  
  ## Reoptimize using same process but on a smaller set of bandwidth values to get more precise
  hseq = seq(hseq[ix-1],hseq[ix+1],length.out = 50)
  sumRes = sapply(hseq,function(h) loocv(dataX,dataY,h,ker))
  
  ## Plot true data with predicted values overlaid for optimized bandwidth, return bandwidth value
  invisible(localLinear(dataX,dataY,hseq[match(min(sumRes),sumRes)],ker))
  return(hseq[match(min(sumRes),sumRes)])
}

loocv <- function(dataX,dataY,h,ker = gaussianKernel){
  #Input:
  # @dataX: a vector of real input data values
  # @dataY: a vector of real output data values
  # @h: bandwith value, single nonnegative value
  # @kernelFunc: kernel function used to evaluate the smoothing weights, defaulted to standard normal
  #Output:
  # @Returns the LOOCV value evaluated for data given, a single nonnegative value that grows with n
  
  ## Find predicted values of y
  yhat = localLinear(dataX,dataY,h,ker)
  
  ## Find weight matrix using wix function, calculate LOOCV errors, sum and return
  weights = sapply(dataX,function(t) wix(dataX,t,h,ker))
  errors <- sapply(1:length(dataX), function(t)  ((dataY[t]-yhat[t])/(1-diag(weights)[t]))^2  )
  return(sum(errors))
}

## Find that perfect bandwidth
hPerf = minH(X,Y,gaussianKernel)

## Heteroscedasticity Check
plot(X,Y-yhat,main='Residual Plot',ylab='Residual Value')

## Since homoscedasticity is not met, we should try to fit our data differently
## Calculate normal predicted y values, but also try using log(Y) and exponentiate results
yhat = localLinear(X,Y,hPerf)
yhatLog = exp(localLinear(X,log(Y),hPerf))

## Plot true values, predicted values from normal smoothing, and predicted values using log scale for Y
## Black = True, Red = Normal, Blue = Log
plot(X,Y)
points(X,yhatLog,col='blue',cex=.5)
lines(X[indices],yhatLog[indices],col='blue')
points(X,yhat,col='red',cex=.5)
lines(X[indices],yhat[indices],col='red')

## Plot log residuals, look more homoscedastic than original
plot(X,log(Y)-log(yhatLog))

## Confidence Bands
## Calculate confidence bounds using fact that fhat is normally distributed
## Calculate weight matrix to help find standard error
w = sapply(X,function(t) wix(X,t,hPerf))
w2 = colSums(w^2)

## Estimate sigma^2
sig2 = sum((Y - yhat)^2/(length(X) - 2*sum(diag(w)) + 2*sum(diag(t(w)%*%w))))

## Find bands and plot along with true data and estimated data
upper = yhat + 1.96 * w2 * sqrt(sig2)
lower = yhat - 1.96 * w2 * sqrt(sig2)
lines(X[indices],upper[indices])
lines(X[indices],lower[indices])

## Try again for log(Y)
yhatLog = localLinear(X,log(Y),hPerf)

upper = yhatLog + 1.96 * w2 * sqrt(sig2)
lower = yhatLog - 1.96 * w2 * sqrt(sig2)

eyhatLog = exp(yhatLog)
eUpper = exp(upper)
eLower = exp(lower)
plot(X,Y)
points(X,eyhatLog,col='blue',cex=.5)
lines(X[indices],eyhatLog[indices],col='red')
lines(X[indices],eUpper[indices])
lines(X[indices],eLower[indices])
