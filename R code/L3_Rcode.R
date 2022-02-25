## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, prompt = TRUE)


## -----------------------------------------------------------------------------
yobs <- c(0,1,1,1,0,1,1,0,1,1)


## -----------------------------------------------------------------------------
# likelihood function: option 1
Lik1.ber <- function(theta, y){
  n = length(y)
  vi = theta^y  * (1-theta)^(1-y)
  return(prod(vi))
}

# log-likelihood function: option 1
lLik1.ber <- function(theta, y){
  n = length(y)
  vi = y * log(theta)  + (1-y)*log(1-theta)
  return(sum(vi))
}


# likelihood function: option 2
Lik2.ber <- function(theta, y){
  vi = dbinom(y, size=1, prob = theta)
  return(prod(vi))
}


## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
theta.gr <- seq(0.1,0.99, len=100)
fLik <- sapply(theta.gr, Lik1.ber,  y=yobs)
flLik <- sapply(theta.gr, lLik1.ber,  y=yobs)
plot(theta.gr, fLik, xlim=c(0,1),lwd=2,
     xlab=expression(theta), 
     type="l",ylab="Likelihood")
plot(theta.gr, flLik, xlim=c(0,1),lwd=2,
     xlab=expression(theta), 
     type="l", ylab="log-likelihood")


## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(theta.gr, fLik/max(fLik), xlim=c(0,1),lwd=2,
     xlab=expression(theta), 
     type="l",ylab="Likelihood")
plot(theta.gr, flLik-max(flLik), xlim=c(0,1),lwd=2,
     xlab=expression(theta), 
     type="l", ylab="log-likelihood")


## -----------------------------------------------------------------------------
# note: by default optimize performs numerical minimisation,
# thus we use the option maximum=TRUE
optimize(f = lLik1.ber, interval = c(0.1, 0.99),
         y = yobs,
         maximum = TRUE)


## -----------------------------------------------------------------------------
# if not already installed
#install.packages("numDeriv")

# after installation, load the package
library(numDeriv)

# check first that the first-order cond. is satisfied
grad(lLik1.ber,0.7, y=yobs)

# observed information at theta=0.7
-hessian(lLik1.ber,0.7, y=yobs)

# check that this agrees with the analytic version
3/(1-0.7)^2 + 7/0.7^2 


## -----------------------------------------------------------------------------
yobs=c(5.1, 7.4, 10.9, 21.3, 12.3, 15.4, 
       25.4, 18.2, 17.4, 22.5)


## -----------------------------------------------------------------------------
# likelihood function
Lik.wei <- function(theta, y){
  # using dwibull(x, alpha, beta)
  vi = dweibull(y, shape = theta[1], scale = theta[2])
  return(prod(vi))
}

# log-likelihood function
lLik.wei <- function(theta, y){
  n = length(y)
  vi = dweibull(y, shape = theta[1], scale = theta[2], log = TRUE)
  
  return(sum(vi))
}


## ---- cache=TRUE--------------------------------------------------------------
# set a regular grid for alpha and beta
theta1 <- seq(0.1, 5, len=300)
theta2 <- seq(12,22, len=300)

# build the bivariate grid
theta12 <- expand.grid(theta1,theta2)

# evaluate the log-Lik over the grid
fLik <-apply(theta12,1,Lik.wei,y=yobs)
flLik <- apply(theta12,1,lLik.wei,y=yobs)

# transform the previous objects in matrices
Lzmat <- matrix(fLik, ncol=300, byrow = F)
lLzmat <-matrix(flLik, ncol=300, byrow = F)

par(mfrow = c(1,2))
image(theta1,theta2, -Lzmat, col = heat.colors(30, rev = TRUE),
      xlab=expression(alpha), 
      ylab=expression(lambda), main="Contours of the Lik")
contour(theta1,theta2, Lzmat, nlevels = 30, add = TRUE)
image(theta1,theta2, -lLzmat,col = heat.colors(30, rev = TRUE),
      xlab=expression(alpha), 
      ylab=expression(lambda),main="Contours of the log-Lik")
contour(theta1,theta2, lLzmat, nlevels = 100, add=TRUE)


## -----------------------------------------------------------------------------
# define the function g(alpha)
g.a <- function(a, y){
  n = length(y)
  oo = sum(log(y)) + n/a - n*sum(y^a * log(y))/sum(y^a)
  
  return(oo)
}


## -----------------------------------------------------------------------------
ga <- sapply(theta1, g.a, y=yobs)
plot(theta1, ga, xlim=c(0.1, 5),
     xlab=expression(alpha), ylab=expression(paste("g(",alpha,")")),
     type="l", lwd=2, col=2)
abline(h=0, lwd=2, lty=2)


## -----------------------------------------------------------------------------
theta1. <- theta1
# set the starting point alpha0
theta0 <- 1

# this is a numericall approximation of dg(a)/da
g.a.prime <- function(a,y) {
  dd <- grad(function(x) g.a(x, y), a)
}


# compute the first 7 iterations
theta1 <- theta0 - g.a(theta0, yobs)/g.a.prime(theta0, yobs)
theta2 <- theta1 - g.a(theta1, yobs)/g.a.prime(theta1, yobs)
theta3 <- theta2 - g.a(theta2, yobs)/g.a.prime(theta2, yobs)
theta4 <- theta3 - g.a(theta3, yobs)/g.a.prime(theta3, yobs)
theta5 <- theta4 - g.a(theta4, yobs)/g.a.prime(theta4, yobs)
theta6 <- theta5 - g.a(theta5, yobs)/g.a.prime(theta5, yobs)
theta7 <- theta6 - g.a(theta6, yobs)/g.a.prime(theta6, yobs)

plot(theta1., ga, xlim=c(0.1, 5),
     xlab=expression(alpha), ylab=expression(paste("g(",alpha,")")),
     type="l", lwd=2, col=2)
abline(h=0, lwd=2, lty=2)

points(y = rep(0,8), 
       x = c(theta0,theta1, theta2,theta3,theta4,theta5,theta6,theta7),
       xlab="iteration", 
       ylab = "hat.alpha at iteration", pch=20, cex=2, 
       type="b", main="Newton-Raphson iterations")
abline(v = 2.75717, lty=2, lwd=2, col="grey")


## -----------------------------------------------------------------------------
(oo <- uniroot(g.a, interval = c(0.1, 10), y=yobs))


## -----------------------------------------------------------------------------
hat.alpha <- oo$root
hat.lambda <- (mean(yobs^hat.alpha))^(1/hat.alpha)


## -----------------------------------------------------------------------------
# we let the legthy output of optim be placed 
# in the object oo and not printed on the screen
# par = c(2,16): is the starting value of the optimizer
# fn = function(x) -lLik.wei(x,yobs) redefines the log-likelihood function 
#  multiplying it by -1.
oo <- optim(par = c(2,1),
            fn = function(x) -lLik.wei(x,yobs),
            method="L-BFGS-B",
            lower = c(0.5,0.1),
            upper = c(4,100),
            hessian = TRUE) # with this option we get the J_n matrix
oo


## -----------------------------------------------------------------------------
eigen(oo$hessian)


## -----------------------------------------------------------------------------
plot(function(x) dweibull(x, shape = oo$par[1], 
                          scale = oo$par[2]), 
     xlim=c(0.01, 45), lwd=2, 
     ylab="Density", xlab="waiting time",
     main= "Estimated probability distribution for the waiting times")
points(x=yobs,y=rep(0,length(yobs)), pch="x")


## -----------------------------------------------------------------------------
oo$par[2]*gamma(1+1/oo$par[1])


## -----------------------------------------------------------------------------
mean(yobs)


## -----------------------------------------------------------------------------
1-pweibull(30, shape = oo$par[1], scale = oo$par[2])


## -----------------------------------------------------------------------------
pweibull(15, shape = oo$par[1], scale = oo$par[2])-
  pweibull(5, shape = oo$par[1], scale = oo$par[2])
