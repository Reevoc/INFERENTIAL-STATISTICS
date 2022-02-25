## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, prompt = TRUE)


## -----------------------------------------------------------------------------------------------------------------------------
bacteria <- read.table("bacteria.csv", sep=";", dec=",", header = T) 

# quick look at the newly created object
str(bacteria)


## -----------------------------------------------------------------------------------------------------------------------------
bacteria$hours <- bacteria$seconds/60^2
with(bacteria, 
     plot(y=A12, hours, ylab="Number of bacteria",xlab = "Time in hours"))


## -----------------------------------------------------------------------------------------------------------------------------
# y = (y1....,y_n)
Y <- as.matrix(bacteria$A12)

# X = [1_n | T | T^2]
X <- as.matrix(cbind(1, bacteria$hours, bacteria$hours^2))

# LS estimate
hat.theta.LS <- solve(t(X) %*% X) %*% t(X) %*% Y
row.names(hat.theta.LS) <- c("hat.theta0", "hat.theta1", "hat.theta2")
hat.theta.LS


## -----------------------------------------------------------------------------------------------------------------------------
with(bacteria,{
  plot(y=A12, hours, ylab="Number of bacteria",xlab = "Time in hours");
  points(hours, 1664.65 + 40.18*hours -2.10* hours^2, type="l", lwd=2,
         lty=2)
})


## -----------------------------------------------------------------------------------------------------------------------------
res <- bacteria$A12 - (1664.65 + 40.18*bacteria$hours -2.10* bacteria$hours^2)
plot(bacteria$hours,res,type="h", ylab="Residuals (Observed - predicted)")


## -----------------------------------------------------------------------------------------------------------------------------
ba.lm <- lm(A12~hours + I(hours^2), data=bacteria)
summary(ba.lm)


## -----------------------------------------------------------------------------------------------------------------------------
# hat.sigma
(hat.sigma = sum(res^2)/length(res))
# square root of hat sigma
sqrt(hat.sigma)


## -----------------------------------------------------------------------------------------------------------------------------
# dividing by n-3 we get the unbiased estimate
(R.hat.sigma = sum(res^2)/(length(res)-3))
sqrt(R.hat.sigma)


## ---- fig.height=10, fig.width=10---------------------------------------------------------------------------------------------
par(mfrow=c(2,2))

# sample size n=5
n <- 5
# exact distribution of n*hat.theta
plot(x= 0:20,y = sapply(0:20, function(x) dpois(x,lambda=3*n/2)), 
     lwd=2,col=2, type="h",xlim=c(0,20), ylab=expression(f(n*hat(lambda)[n])), 
     xlab = expression(n*hat(lambda)[n]), main=paste("n =",n))
# approximate dist. of n*hat.theta
plot(function(x) dnorm(x, mean=n*3/2, sd=sqrt(3*n/2)), 
     lwd=2, add=TRUE, xlim = c(0,20))
legend("topright", legend = c("exact", "asymptotic"), lwd=c(2,2), col=c(2,1),
       lty=c(1,1))

n <- 10
plot(x= 0:35,y = sapply(0:35, function(x) dpois(x,lambda=3*n/2)), 
     lwd=2,col=2, type="h",xlim=c(0,35), ylab=expression(f(n*hat(lambda)[n])), 
     xlab = expression(n*hat(lambda)[n]),main=paste("n =",n))
plot(function(x) dnorm(x, mean=n*3/2, sd=sqrt(3*n/2)), 
     lwd=2, add=TRUE, xlim = c(0,35))


n <- 30
plot(x= 20:80,y = sapply(20:80, function(x) dpois(x,lambda=3*n/2)), 
     lwd=2,col=2, type="h",xlim=c(20,80), ylab=expression(f(n*hat(lambda)[n])), 
     xlab = expression(n*hat(lambda)[n]),main=paste("n =",n))
plot(function(x) dnorm(x, mean=n*3/2, sd=sqrt(3*n/2)), 
     lwd=2, add=TRUE, xlim = c(20,80))

n <- 50
plot(x= 40:110, y = sapply(40:110, function(x) dpois(x,lambda=3*n/2)), 
     lwd=2,col=2, type="h",xlim=c(40,110), ylab=expression(f(n*hat(lambda)[n])), 
     xlab = expression(n*hat(lambda)[n]),main=paste("n =",n))
plot(function(x) dnorm(x, mean=n*3/2, sd=sqrt(3*n/2)), 
     lwd=2, add=TRUE, xlim = c(40,120))


## -----------------------------------------------------------------------------------------------------------------------------
# sample of size 5
n1 = 5
plot(function(x) dnorm(x, mean=3/2, sd=sqrt(3/(2*n1))), 
     lwd=2, xlim=c(0,3.5),n=500, ylim=c(0,2.5),
     ylab=expression(f(hat(lambda)[n])), 
     xlab = expression(hat(lambda)[n]),
     main = "Distribution of the MLE")
abline(v=3/2, lwd=2, lty=2, col="lightgray")
# samples of sizes 10, 30, 50, 
sapply(c(10, 30, 50), 
       function(y) plot(function(x) dnorm(x, mean=3/2, sd=sqrt(3/(2*y))), 
                        lwd=2, xlim=c(0,10),col=round(1+y/10),lty=round(1+y/10), n=500, add=TRUE))
legend("topright", legend = c("n=5", "n=10", "n=30", "n=50"),
       lwd=c(2,2,2,2), col=c(1,2,4,6), lty=c(1,2,4,6))


## -----------------------------------------------------------------------------------------------------------------------------
challenger <- read.table("challenger2.dat", header = TRUE)
plot(challenger)
points(x=53, y=5, col=2, pch="*", cex=2)


## -----------------------------------------------------------------------------------------------------------------------------
plot(challenger[challenger$orings>0,],
     ylab = "Orings(>0)")
points(x=53, y=5, col=2, pch="*", cex=2)


## -----------------------------------------------------------------------------------------------------------------------------
# log-likelihood function for a logistic regression with a single covariate
logLBin <- function(theta, data){
  
  n  = nrow(data)
  alpha. = theta[1]
  beta. = theta[2] 
  
  # compute the success prob. for each launch (observation), i.e. e^(a + b*t)/(1+e^(a+b*t))
  thetai = plogis(alpha. + beta.*data$temperature)
  
  # log-likelihood for each launch (observation)
  oo = dbinom(x = data$orings, # observed broken o-rings
              size = 6, # maximum number of o-rings
              prob = thetai, 
              log=TRUE # TRUE gives the log of the density
  )
  
  # sum all the log-likelihood contributions
  return(sum(oo))
}


## -----------------------------------------------------------------------------------------------------------------------------
# actually optim minimises the negative log-likelihood
start0 <- c(10, 0)
(oo <- optim(par = start0, 
             fn = function(x) -logLBin(x, data = challenger), 
             hessian = TRUE)
)


## -----------------------------------------------------------------------------------------------------------------------------
# estimated prob. of failure of single O-ring
(theta.t <- plogis(oo$par[1] + oo$par[2]*31))

# estimated prob. distribution of nr. of failed O-rings at temp = 31F
plot(x = 0:6, dbinom(0:6, 6, theta.t), type="h", 
     xlab="O-rings", ylab="p.d.f.")


## -----------------------------------------------------------------------------------------------------------------------------
sum(dbinom(4:6,6,theta.t))

# or also
1-pbinom(3,6, prob = theta.t,)
