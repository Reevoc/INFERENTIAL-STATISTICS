## ----setup, include=FALSE--------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, prompt = TRUE)


## --------------------------------------------------------------------------------------------
set.seed(5) # fix the random seed
n <- 10 # set the sample size
mu0 <-  0
sigma2.0 <- 2
yobs <- rnorm(n, mu0, sqrt(sigma2.0))
yobs # measurements we are given by Nature


## --------------------------------------------------------------------------------------------
alpha = 0.05
(bar.y <- mean(yobs)) # sample mean
(se <- sqrt(sigma2.0/n)) # this is the standard error of the estimate
bar.y + c(-1,1)*qnorm(p = alpha/2, lower.tail =FALSE)*se
# note: we used lower.tail =FALSE in the quantile function of the
# normal distribution in order to have the upper quantile as desire.



## ---- cache=TRUE-----------------------------------------------------------------------------
set.seed(5)
N = 1e+6 # the number of intervals we want to calculate
my.N.CI <- matrix(NA, nrow = N, ncol = 2)# put the obs.intervals here
for(i in 1:N){
  yi <- rnorm(n, mu0, sqrt(sigma2.0))
  bar.yi <- mean(yi)
  my.N.CI[i,] <- bar.yi + c(-1,1)*qnorm(p = alpha/2, lower.tail =FALSE)*sqrt(sigma2.0/n)
}
# here are the first 6 observed intervals out of N
head(my.N.CI)


## --------------------------------------------------------------------------------------------
# we use the plotrix library for doing this plot
# if not already installed use
# install.packages(plotrix) 
library(plotrix)
plotCI(x= 1:50, y = apply(my.N.CI[1:50,], 1, mean), 
       li = my.N.CI[1:50,1],ui = my.N.CI[1:50,2],
       xlab="Observed confidence intervals for mu",ylab=NA, lwd=3,)
abline(h = 0, lwd=2, lty=2)


## --------------------------------------------------------------------------------------------
mu0.inside <- apply(my.N.CI, MARGIN = 1,
                    function(x) ifelse(mu0 >= x[1]  & mu0 <= x[2],1,0))
# m0.inside is a vector of 1 and 0
head(mu0.inside)

# how many ones are there relative to N?
mean(mu0.inside)


## --------------------------------------------------------------------------------------------
set.seed(5) # fix the random seed
n <- 10 # set the sample size
mu0 <-  0
sigma2.0 <- 1
yobs <- rnorm(n, mu0, sqrt(sigma2.0))


## --------------------------------------------------------------------------------------------
(hat.sig2.mu <- sum((yobs-mu0)^2)/n)
# CI is
c(n*hat.sig2.mu/qchisq(p=alpha/2, df=n, lower.tail = FALSE),
  n*hat.sig2.mu/qchisq(p=1-alpha/2, df=n, lower.tail = FALSE))


## --------------------------------------------------------------------------------------------
set.seed(5) # fix the random seed
n <- 10 # set the sample size
mu0 <-  0
sigma2.0 <- 1
yobs <- rnorm(n, mu0, sqrt(sigma2.0))


## --------------------------------------------------------------------------------------------
# recall that var divides by n-1!
(hat.sig2 <- var(yobs))
(hat.mu <- mean(yobs))

# the standard error of hat.mu
se <- sqrt(hat.sig2/n)

# CI for mu
c(hat.mu + c(-1,1)*qt(alpha/2, df=n-1, low=F)*se)

# CI for mu
c(n*hat.sig2/qchisq(p=alpha/2, df=n-1, lower.tail = FALSE),
  n*hat.sig2/qchisq(p=1-alpha/2, df=n-1, lower.tail = FALSE))


## --------------------------------------------------------------------------------------------
# for the CI for mu
t.test(yobs)


## --------------------------------------------------------------------------------------------
# read the file, specifying header=TRUE since the first row contains the names of the variables.
motors <- read.table("wm_motors.txt", header = TRUE)
head(motors)

# another quick view for data frames is
str(motors)


## --------------------------------------------------------------------------------------------
# energy consumption with motor GEN1
y <- motors$energy[motors$motor == "GEN1"]

# energy consumption with motor GEN2
x <- motors$energy[motors$motor == "GEN2"]


## --------------------------------------------------------------------------------------------
# sample averages of energy consumption under GEN1 and GEN2
bar.y <- mean(y)
bar.x <- mean(x)

# compute the size of the two samples
n <- length(y)
m <- length(x)

# sample variances 
s2.y <- var(y)
s2.x <- var(x)

# pooled variance
pooled.s2 <- ((n-1)*s2.x + (m-1)*s2.y)/(n+m-2)

# the standard error is thus
se = sqrt(pooled.s2*(1/n+1/m))

# the confidence interval is then
(bar.y-bar.x) + c(-1,1)*qt(alpha/2, df=n+m-2, low=F)*se


## --------------------------------------------------------------------------------------------
t.test(y,x, var.equal = T)


## --------------------------------------------------------------------------------------------
c(s2.y/s2.x * qf(alpha/2, df1=n-1, df2=m-1),
  s2.y/s2.x * qf(1- alpha/2, df1=n-1, df2=m-1))

# the same done with the built-in R command
var.test(y,x)


## --------------------------------------------------------------------------------------------
lambda0 <- 3/2
N <- 1e+6

# we first build a function which takes as input an observed sample and
# outputs a confidence interval
ci.poisson <- function(ysamp, alpha){
  n = length(ysamp)
  bar.y = mean(ysamp)
  se = sqrt(bar.y/n) 
  oo = c(bar.y + c(-1,1)*qnorm(alpha/2, low=F)*se)
  return(oo)
}


## ---- cache=TRUE-----------------------------------------------------------------------------
set.seed(5)
n <- 5
CI.5 <- matrix(NA, nrow = N, ncol = 2)# put the obs.intervals here
for(i in 1:N) {
  y5 <- rpois(n, lambda0)
  CI.5[i,] <- ci.poisson(y5,alpha = 0.05)
}
lambda0.inside5 <- apply(CI.5, MARGIN = 1,
                         function(x) ifelse(lambda0 >= x[1]  & lambda0 <= x[2],1,0))
# how many ones are there?
mean(lambda0.inside5)


## ---- cache=TRUE-----------------------------------------------------------------------------
set.seed(5)
n <- 50
CI.50 <- matrix(NA, nrow = N, ncol = 2)# put the obs.intervals here
for(i in 1:N) {
  y50 <- rpois(n, lambda0)
  CI.50[i,] <- ci.poisson(y50,alpha = 0.05)
}
lambda0.inside50 <- apply(CI.50, MARGIN = 1,
                          function(x) ifelse(lambda0 >= x[1]  & lambda0 <= x[2],1,0))
# how many ones are there?
mean(lambda0.inside50)


## --------------------------------------------------------------------------------------------
# 36 green out of 56
bar.y = 39/(39+17)
se = sqrt((bar.y)*(1-bar.y)/56)

# the 0.95 CI is
bar.y + c(-1,1)*qnorm(alpha/2, low=F)*se


## --------------------------------------------------------------------------------------------
prop.test(x=39, n=39+17, correct = F)


## --------------------------------------------------------------------------------------------
challenger <- read.table("challenger2.dat", header = TRUE)
plot(challenger)
points(x=53, y=5, col=2, pch="*", cex=2)


## --------------------------------------------------------------------------------------------
plot(challenger[challenger$orings>0,],
     ylab = "Orings(>0)")
points(x=53, y=5, col=2, pch="*", cex=2)


## --------------------------------------------------------------------------------------------
logLBin <- function(theta, data){
  
  n  = nrow(data)
  alpha. = theta[1]
  beta. = theta[2] 
  
  # compute the success prob. for each launch (observation)
  thetai = plogis(alpha. + beta.*data$temperature)
  
  # log-likelihood for each launch (observation)
  oo = dbinom(x = data$orings, size = 6, prob = thetai, log=TRUE)
  
  # sum all the log-likelhood contributions
  return(sum(oo))
}


## --------------------------------------------------------------------------------------------
start0 <- c(10, 0)
(oo <- optim(par = start0, fn = function(x) -logLBin(x, data = challenger), 
             hessian = TRUE))


## --------------------------------------------------------------------------------------------
var.theta <- solve(oo$hessian)
(se.beta <- sqrt(var.theta[2,2]))
oo$par[2] + c(-1,1)*qnorm(alpha/2, low=F)*se.beta


## --------------------------------------------------------------------------------------------
# estimated probability of "success"
(theta.t <- plogis(oo$par[1] + oo$par[2]*31))
plot(x = 0:6, dbinom(0:6, 6, theta.t), type="h", 
     xlab="O-rings", ylab="p.d.f.")


## --------------------------------------------------------------------------------------------
sum(dbinom(4:6,6,theta.t))

# or also by
1-pbinom(3,6, prob = theta.t,)


## --------------------------------------------------------------------------------------------
oo2 <- glm(cbind(orings, 6-orings)~temperature, data=challenger,
           family = binomial())
summary(oo2)
