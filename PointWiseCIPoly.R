#Pointwise Confidence Interval Plot

PWConfPlot <- function(xData,yDataMat,confLevel,colorReq='red'){
  numPoints = length(xData)
  low = (1 - confLevel)/2
  high = 1 - low
  bounds = sapply(1:numPoints,function(t) quantile(yDataMat[t,],c(low,high)))
  avgVals = sapply(1:numPoints,function(t) mean(yDataMat[t,]))
  rx = rev(xData)
  rb = rev(bounds[2,])
  x = c(xData,rx)
  y = c(bounds[1,],rb)
  #plot(xData,avgVals, col='red',type = 'l',lty='dashed')
  polygon(x,y,density=100,col=colorReq)
  #lines(xData,avgVals,col='red',lty='dashed',lwd=2)
}

#PWConfPlot(data,simData,.95)