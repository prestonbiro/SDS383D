#SDS 386
#Exercise 3
#Kernal Smoothing
setwd("~/Documents/SDS 383D")
data = read.csv('Data/gdpgrowth.csv')

weightFunc <- function(prevX,newX,kerFunc,bandw = .1){
  #Input:
  #prevX is a vector of previous data points x
  #newX is a real valued individual new x point which we wish to estimate a corresponding
  #   y point (not calculated here, but weights calculated to take a step towards it)
  #kerFunc is a smooth function satisfying three properties (E[1] = 1, E[x] = 0,
  #   E[x^2] >0 all with respect to kerFunc)
  #bandw is the bandwidth being used, default to .1
  #Output:
  #weightVec is a vector of weights corresponding to each of the data points prevX
  h = bandw
  n = length(prevX)
  weightVec = rep(NA,n)
  for(i in 1:n){
    dist = prevX[i] - newX
    weightVec[i] = 1/h * kerFunc(dist,h)
  }
  return(weightVec)
}

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

uniformKernel <- function(distance,bandwidth){
  #Input:
  #distance is a real number (can be negative, but will take magnitude)
  #bandwidth is a positive real number 
  #Output:
  #The distance divided by the bandwidth evaluated by the Uniform kernel function
  #   K(x) = 1/2 if abs(x) <= 1, 0 else
  #eval is the nonnegative real value returned for the evaluated value
  dist = abs(distance)
  xEval = dist/bandwidth
  if(xEval <= 1) eval = 1/2
  else eval = 0
  return(eval)
}

linearSmoother <- function(prevX,newX,realY,kernelFunc,weightFinder,bandw = .1){
  #Input:
  #prevX is a vector of previous data points x
  #newX is a real valued individual new x point which we wish to estimate a corresponding
  #   y point (not calculated here, but weights calculated to take a step towards it)
  #realY is the vector of real output values corresponding to prevX
  #kerFunc is a smooth function satisfying three properties (E[1] = 1, E[x] = 0,
  #   E[x^2] >0 all with respect to kerFunc)
  #weightFinder is a function that takes inputs of the prevX, newX, and kernalFunc to 
  #   find the weight values for the individual prevX values
  #bandw is the bandwidth being used, default to .1
  #linear smoother will send the prevX, newX, and kernelFunc to weightFinder and 
  #   use the weights to make predict the output of a function evaluated at newX using
  #   f_hat(newX) = sum_i[ weights(prevX_i,newX) * realY_i ]
  #Output:
  #fhat is the real valued estimate of the function evaluated at newX
  
  weights = weightFinder(prevX,newX,kernelFunc,bandw)
  normWeights = weights/sum(weights)
  n = length(realY)
  indEst = rep(NA,n)
  for(i in 1:n){
    indEst[i] = normWeights[i] * realY[i]
  }
  fhat = sum(indEst)
  return(fhat)
}

X = data$LGDP60
Y = data$LIFE60
B = 500
xnewSeq = seq(6,9,length.out=B)
yu = rep(NA,B)
yg = rep(NA,B)
for(i in 1:B){
  yu[i] = linearSmoother(X,xnewSeq[i],Y,uniformKernel,weightFunc,1)
  yg[i] = linearSmoother(X,xnewSeq[i],Y,gaussianKernel,weightFunc,1)
}
plot(X,Y)
points(xnewSeq,yu,col='red')
points(xnewSeq,yg,col='blue')


nonLinFunc <- function(x) exp(-x^2) + x^2 -x 
curve(nonLinFunc,xlim=c(-3,3))
n = 100
Xseq = seq(-3,3,length.out = n)
noiseX = rnorm(n)
Xrand = Xseq + noiseX
Xrand = Xrand - mean(Xrand)
noiseY = rnorm(n)
Yrand = nonLinFunc(Xrand) + noiseY
Yrand = Yrand - mean(Yrand)
plot(Xrand,Yrand)

B = 200
xStars = seq(-2.5,2.5,length.out = B)
yUnif = rep(NA,B)
yGaus = rep(NA,B)
for(i in 1:B){
  yUnif[i] = linearSmoother(Xrand,xStars[i],Yrand,uniformKernel,weightFunc,.01)
  yGaus[i] = linearSmoother(Xrand,xStars[i],Yrand,gaussianKernel,weightFunc,.01)
}
points(xStars,yUnif,col='red')
points(xStars,yGaus,col='blue')
#curve(nonLinFunc,add=T)

smoothWrapper <- function(trainX,trainY,testX,testY,smoothFunction = linearSmoother,weightCalculator = weightFunc,kernelFunction = gaussianKernel){
  #Input:
  #trainX = training data set vector of real inputs X
  #trainY = training data set vector of real inputs X
  #testX = testing data set vector of real inputs Y
  #testY = testing data set vector of real inputs Y
  #Other inputs are functions used to 
  hSeq = seq(.01,1,length.out = 100)
  x0 = trainX - mean(trainX)
  y0 = trainY - mean(trainY)
  xt = testX - mean(testX)
  yt = testY - mean(testY)
  yhat = matrix(NA,nrow=length(testX),ncol=length(hSeq))
  for(j in 1:length(hSeq)){
    for(i in 1:length(testX)){
      yhat[i,j] = smoothFunction(x0,xt[i],y0,kernelFunction,weightCalculator,hSeq[j])
    }
  }
  
  eps2 = (yhat - yt)^2
  hIdeal = hSeq[match(min(colSums(eps2)),colSums(eps2))]
  print(hIdeal)
  hgood = hSeq[(match(hIdeal,hSeq)-1):(match(hIdeal,hSeq)+2)]
  par(mfrow=c(2,2))
  for(j in 1:4){
    plot(x0,y0,main= hgood[j])
    #points(xt,yt,pch=2,col=5*j)
    points(xt,yhat[,j],pch=8,col=3*j)
  }
  
  return(yhat[,match(min(colSums(eps2)),colSums(eps2))])
}


X = data$LGDP60
Y = data$LIFE60
trainingX = X[1:ceiling(.75 * length(X))]
testingX = X[(ceiling(.75 * length(X))+1):length(X)]
trainingY = Y[1:ceiling(.75 * length(Y))]
testingY = Y[(ceiling(.75 * length(Y))+1):length(Y)]
smoothWrapper(trainingX,trainingY,testingX,testingY)


#2 by 2 Table Scenario
par(mfrow=c(1,1))
N = 500
wigglyFunc <- function(x) 4 * cos(30 * x)
curve(wigglyFunc)
smoothFunc <- function(x) 4 - 16*(x-.5)^2
curve(smoothFunc)
x1 = runif(N)
x2 = runif(N)
x3 = runif(N)
x4 = runif(N)
par(mfrow=c(2,2))
quietWiggles = .5*rnorm(N) + wigglyFunc(x1)
noisyWiggles = 1.5*rnorm(N) + wigglyFunc(x2)
quietSmooth = .25*rnorm(N) + smoothFunc(x3)
noisySmooth = 1.25*rnorm(N) + smoothFunc(x4)
plot(x1,quietWiggles)
plot(x2,noisyWiggles)
plot(x3,quietSmooth)
plot(x4,noisySmooth)

train1x = x1[1:(.8*N)]
train2x = x2[1:(.8*N)]
train3x = x3[1:(.8*N)]
train4x = x4[1:(.8*N)]
test1x = x1[(.8*N + 1):N]
test2x = x2[(.8*N + 1):N]
test3x = x3[(.8*N + 1):N]
test4x = x4[(.8*N + 1):N]
train1y = quietWiggles[1:(.8*N)]
train2y = noisyWiggles[1:(.8*N)]
train3y = quietSmooth[1:(.8*N)]
train4y = noisySmooth[1:(.8*N)]
test1y = quietWiggles[(.8*N + 1):N]
test2y = noisyWiggles[(.8*N + 1):N]
test3y = quietSmooth[(.8*N + 1):N]
test4y = noisySmooth[(.8*N + 1):N]

yhat1 = smoothWrapper(train1x,train1y,test1x,test1y)
yhat2 = smoothWrapper(train2x,train2y,test2x,test2y)
yhat3 = smoothWrapper(train3x,train3y,test3x,test3y)
yhat4 = smoothWrapper(train4x,train4y,test4x,test4y)
par(mfrow=c(2,4))
plot(train1x,train1y)
plot(train2x,train2y)
plot(train3x,train3y)
plot(train4x,train4y)
plot(test1x,yhat1)
plot(test2x,yhat2)
plot(test3x,yhat3)
plot(test4x,yhat4)















