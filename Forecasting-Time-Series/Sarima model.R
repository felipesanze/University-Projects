library(fBasics)
library(forecast) 

# Set the folder (path) that contains this R file as the working directory
dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)

# *** SESSION 8. Real Example ***
datos<-read.csv("Data.csv",header=TRUE,sep=";",dec=",")
y<-datos[,2] 

ts.plot(y)

par(mar=c(5,5,5,5))

nlags=80   # play with this parameter..
par(mfrow=c(2,1))
acf(y,nlags)
pacf(y,nlags)  

s=4       # seasonal parameter FOR THIS DATA SET

ndiffs(y, alpha=0.05, test=c("adf")) # regular differences?

nsdiffs(y,m=s,test=c("ocsb"))  # seasonal differences?

z<-diff(y)

ts.plot(z)  
par(mfrow=c(2,1))
acf(z, nlags)  
pacf(z, nlags)

ndiffs(z, alpha=0.05, test=c("adf"))

# estimate the SAR and analyze the estimated parameters. Compare with the Seasonal Difference
fit<-arima(y,order=c(1,1,0),seasonal=list(order=c(1,0,0),period=s)) 
fit

ts.plot(fit$residuals)
par(mfrow=c(2,1))
acf(fit$residuals,nlags)
pacf(fit$residuals,nlags)    

ndiffs(fit$residuals, alpha=0.05, test=c("adf")) # regular differences?

Box.test(fit$residuals,lag=25) # white noise residuals?

Box.test(fit$residuals,lag=48)   # play with the number of lags # a demander à raquel

shapiro.test(fit$residuals)  # normality test


# normality test graphically
par(mfrow=c(1,1))
hist(fit$residuals,prob=T,ylim=c(0,1),xlim=c(mean(fit$residuals)-3*sd(fit$residuals),mean(fit$residuals)+3*sd(fit$residuals)),col="red")
lines(density(fit$residuals),lwd=2)
mu<-mean(fit$residuals)
sigma<-sd(fit$residuals)
x<-seq(mu-3*sigma,mu+3*sigma,length=100)
yy<-dnorm(x,mu,sigma)
lines(x,yy,lwd=2,col="blue")


# squared residuals

ts.plot(fit$residuals^2)
par(mfrow=c(2,1))
acf(fit$residuals^2,nlags)
pacf(fit$residuals^2,nlags)    

Box.test(fit$residuals^2,lag=25) # SWN

# point predictions and standard errors

y.pred<-predict(fit,n.ahead=24)
y.pred$pred   # point predictions
y.pred$se    # standard errors

shapiro.test(y.pred$pred)

# plotting real data with point predictions

new <- c(y,y.pred$pred) # real data + predicted values

plot.ts(new,main="Predictions",
        ylab="Dollars",col=3,lwd=2) # time series plot
lines(y,col=4,lwd=2) # for the second series
legend("topleft",legend=c("Predictions","Historical"),col=c(3,4),
       bty="n",lwd=2)
