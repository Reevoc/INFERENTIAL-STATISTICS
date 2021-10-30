## ----setup, include=FALSE------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, prompt = TRUE)


## ------------------------------------------------------------------------------------------------------------
head(cars)
summary(cars)


## ------------------------------------------------------------------------------------------------------------
var(cars$speed)
var(cars$dist)


## ------------------------------------------------------------------------------------------------------------
median(cars$speed)
median(cars$dist)


## ------------------------------------------------------------------------------------------------------------
quantile(cars$speed, probs = 1/4, type=6)# lower quarile
quantile(cars$speed, probs = 3/4, type=6) # upper quartile
quantile(cars$speed, probs = 2/4, type=6) # the median


## ------------------------------------------------------------------------------------------------------------
x <- c(0.318, 1.765, 0.259, 0.450, 0.730, 0.235, 0.017, 1.010, 1.418, 0.480)
(x.order <- sort(x))
n <- length(x)
k <- (1:n)/(n+1)
qq <- -2*log(1-k)/3
plot(x=x.order, y=qq, xlab="Observed quantiles", ylab="Theoretical quantiles")
abline(a=0, b=1)


## ------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
qqnorm(cars$speed,main="QQ-plot: normal vs speed")
qqline(cars$speed)

qqnorm(cars$dist,main="QQ-plot: normal vs dist")
qqline(cars$dist)


## ------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
# a non scaled histogram
hist(cars$speed)
# a scaled histogram
hist(cars$speed, freq = FALSE)


## ------------------------------------------------------------------------------------------------------------
table(cars$speed)


## ------------------------------------------------------------------------------------------------------------
plot(0,xlim=c(0,30), ylim=c(0,1), cex=0, ylab="Observed Empirical d.f.")
segments(x0=0,x1=4,y0=0,y1=0, lwd=2)

segments(x0=4,x1=7,y0=2/50,y1=2/50, lwd=2)
points(x=4, y=2/50, pch=20, cex=2)

segments(x0=7,x1=8,y0=4/50,y1=4/50, lwd=2)
points(x=7, y=4/50, pch=20, cex=2)

# etc
# etc

segments(x0=25,x1=30,y0=50/50,y1=50/50, lwd=2)
points(x=25, y=50/50, pch=20, cex=2)


## ------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(ecdf(cars$speed))
plot(ecdf(cars$dist))


## ------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))
boxplot(cars$speed, main="Boxplot of speed")
boxplot(cars$dist, main="Boxplot of dist")


## ------------------------------------------------------------------------------------------------------------
set.seed(2020)
m = 1e+5
n <- 10
# to be filled by exponential r.variates
x.exp <- vector(mode = "numeric", n)

# to be filled by draws from X_{(n)}
rdraw <- vector(mode = "numeric", m)

for(i in 1:m) {
  x.exp <- rexp(n,rate=1)
  
  # sort exponential r.variates
  x.exp.order <- sort(x.exp)
  
  # save the largest 
  rdraw[i] <- x.exp.order[n]
}


## ------------------------------------------------------------------------------------------------------------
ff <- function(x,n) n*(1-exp(-x))^(n-1)*exp(-x)

hist(rdraw, freq = F, breaks = 40)
plot(function(x) ff(x,n=n), xlim=c(0,10), n=100, add=T, lwd=2)


## ------------------------------------------------------------------------------------------------------------
plot(cars$speed, cars$dist)


## ------------------------------------------------------------------------------------------------------------
# covariance
cov(cars$speed, cars$dist)

# correlation: option 1
cor(cars$speed, cars$dist)

# correlation: option 2
cov(cars$speed, cars$dist)/
  sqrt(var(cars$speed)*var(cars$dist))
