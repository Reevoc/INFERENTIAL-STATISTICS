## ----setup, include=FALSE---------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, prompt = TRUE)


## ---------------------------------------------------------------------------------------
a <- 1:10


## ---------------------------------------------------------------------------------------
a[3]


## ---------------------------------------------------------------------------------------
a[1:3]


## ---------------------------------------------------------------------------------------
a[c(2,3,7)]


## ---------------------------------------------------------------------------------------
a <- c(1,2,3,4,5,6,7,8,9,10)


## ---------------------------------------------------------------------------------------
A <- matrix(c(1:10), ncol=2)


## ---------------------------------------------------------------------------------------
# (1,2) element of A is
A[1,2]


## ----pdf-unif, fig.cap="The p.d.f (left) and d.f. of the Unif(0,1) r.v.."---------------
par(mfrow = c(1,2))
plot(dunif, xlim=c(-1,2),
     ylab="f(x)", main = "The p.d.f of the Unif(0,1) r.v.")
plot(punif, xlim=c(-1,2),
     ylab="F(x)", main = "The d.f of the Unif(0,1) r.v.")


## ---------------------------------------------------------------------------------------
set.seed(2020)
xunif <- runif(1e+4)


## ---------------------------------------------------------------------------------------
par(mfrow = c(1,1))
hist(xunif, xlim=c(-1,2), freq= F)
plot(dunif, xlim=c(-1,2), add=TRUE, col=2, lwd=2)


## ---------------------------------------------------------------------------------------
set.seed(2020)
xunif <- runif(3e+5)
hist(xunif, xlim=c(-1,2), freq = FALSE)
plot(dunif, xlim=c(-1,2), add=TRUE, col=2, lwd=2)


## ---------------------------------------------------------------------------------------
set.seed(2020)
xunif1 <- runif(10)
set.seed(2020)
xunif2 <- runif(50)
set.seed(2020)
xunif3 <- runif(1e+2)
set.seed(2020)
xunif4 <- runif(5e+2)
set.seed(2020)
xunif5 <- runif(1e+3)
averages <- c(mean(xunif1), mean(xunif2), mean(xunif3),mean(xunif4),mean(xunif5))
plot(x = c(10,50,1e+2,5e+2,1e+3), y=averages, ylim=c(0.3, 0.6), pch=20,
     ylab = "Sample average",
     xlab="sample size")
abline(h = 1/2, lwd=2, lty=2)


## ---------------------------------------------------------------------------------------
set.seed(2020)
# generate 1000 r. draws from the N(0,1) dist
xnorm <- rnorm(1e+3)


## ---------------------------------------------------------------------------------------
# generate 1000 r. draws from the N(1/2,2) dist
xnorm1 <- rnorm(1e+3, mean = 1/2, sd = sqrt(2))


## ---------------------------------------------------------------------------------------
mu <- 1/2
st.dev <- sqrt(2)
xnorm2 <- mu + st.dev*xnorm


## ---------------------------------------------------------------------------------------
par(mfrow = c(1,2))
hist(xnorm1, xlim=c(-5,5), freq = FALSE)
plot(function(x) dnorm(x, mean = 1/2, sd=sqrt(2)),
     xlim=c(-5,5), add=TRUE, col=2, lwd=2)
hist(xnorm2, xlim=c(-5,5), freq = FALSE)
plot(function(x) dnorm(x, mean = 1/2, sd=sqrt(2)),
     xlim=c(-5,5), add=TRUE, col=2, lwd=2)



## ---------------------------------------------------------------------------------------
fx <- function(x) {
  out1 <- x^3 + log(x)
  out2 <- out1/x
  return(out2)
}


## ---------------------------------------------------------------------------------------
set.seed(2020)
znorm <- rnorm(1e+4)
ychi <- znorm^2
hist(ychi, freq = FALSE)
plot(function(x) dchisq(x, df=1), add=TRUE, lwd=2, col=2, xlim=c(0,20))


## ---------------------------------------------------------------------------------------
set.seed(2020)
znorm1 <- rnorm(1e+4)
znorm2 <- rnorm(1e+4)
ychi2 <- znorm1^2+ znorm2^2
hist(ychi2, freq = FALSE)
plot(function(x) dchisq(x, df=2), add=TRUE, lwd=2, col=2, xlim=c(0,20))


## ---------------------------------------------------------------------------------------
# drawing from t-Student distribution with nu=4
# using its definition
set.seed(2020)
nu <- 4
znorm <- rnorm(1e+4)
xchi <- rchisq(1e+4, df=4)
yt <- znorm/sqrt(xchi/nu)


## ---------------------------------------------------------------------------------------
hist(yt, freq = FALSE, breaks = 50, ylim=c(0,0.45))
plot(function(x) dt(x, df=nu), add=TRUE, lwd=2, col=2, xlim=c(-10,10))


## ---------------------------------------------------------------------------------------
set.seed(2020)
unif <- runif(1e+4)
yexp <- -log(1-unif)
# alternatively we can also use
# yexp1 <- -log(unif)
# becasue U ~ Unif(0,1), then also 1-U ~ Unif(0,1)


## ---------------------------------------------------------------------------------------
hist(yexp, freq = FALSE, breaks = 50)
plot(dexp, add=TRUE, lwd=2, col=2, xlim=c(0,10))


## ---------------------------------------------------------------------------------------
par(mfrow=c(2,2))
plot(x = 0:10, dbinom(0:10, size=10, prob=1/2), type="h",
     ylab="f(x)", xlab="x", main="p.d.f. of the Bin(10,1/2) r.v.")
plot(x = 0:10, dbinom(0:10, size=10, prob=1/10), type="h",
     ylab="f(x)", xlab="x", main="p.d.f. of the Bin(10,1/10) r.v.")
plot(x = 0:20, dbinom(0:20, size=20, prob=1/2), type="h",
     ylab="f(x)", xlab="x", main="p.d.f. of the Bin(20,1/2) r.v.")
plot(x = 0:20, dbinom(0:20, size=20, prob=1/10), type="h",
     ylab="f(x)", xlab="x", main="p.d.f. of the Bin(20,1/10) r.v.")


## ---------------------------------------------------------------------------------------
par(mfrow = c(2,2))
plot(function(x) pbinom(x, size=10, prob=1/2),xlim=c(0,10), n=500,
     ylab="F(x)", main = " d.f. of the Bin(10,1/2) r.v.")
plot(function(x) pnorm(x, mean=10*1/2, sd=sqrt(10*1/4)),xlim=c(0,10),
     add=TRUE, col=2, lwd=2, lty=2)

plot(function(x) pbinom(x, size=10, prob=1/10),xlim=c(0,10), n=500,
     ylab="F(x)", main = " d.f. of the Bin(10,1/10) r.v.")
plot(function(x) pnorm(x, mean=10*1/10, sd=sqrt(10*9/100)),xlim=c(0,10),
     add=TRUE, col=2, lwd=2, lty=2)


plot(function(x) pbinom(x, size=20, prob=1/2),xlim=c(0,20), n=500,
     ylab="F(x)", main = " d.f. of the Bin(20,1/2) r.v.")
plot(function(x) pnorm(x, mean=20*1/2, sd=sqrt(20*1/4)),xlim=c(0,20),
     add=TRUE, col=2, lwd=2, lty=2)

plot(function(x) pbinom(x, size=20, prob=1/10),xlim=c(0,20), n=500,
     ylab="F(x)", main = " d.f. of the Bin(20,1/10) r.v.")
plot(function(x) pnorm(x, mean=20*1/10, sd=sqrt(20*9/100)),xlim=c(0,20),
     add=TRUE, col=2, lwd=2, lty=2)


## ---------------------------------------------------------------------------------------
# real-valued function with arguments:
# x: vector (2 x 1), all other argumens are scalars
dbvnorm <- function(x, mu1, mu2, sigma1.2, sigma2.2, sigma12){
  # build the covariance matrix Sigma
  Sigma = matrix(c(sigma1.2, sigma12, sigma12,sigma2.2),ncol=2)
  
  # build the vector of means
  mu = c(mu1, mu2)
  
  # compute a first part the of p.d.f, aka the "kernel"
  out1 = exp(-0.5 * t(x-mu) %*% solve(Sigma) %*% (x-mu))
  
  # scale the kernel by the normalising constant
  out2 = out1/((2*pi)^(2/2) * sqrt(det(Sigma)))
  
  # return the output
  return(out2)
  
  # you can also omit print() by writing just
  # out2
}


## ---------------------------------------------------------------------------------------
mu1 <- 1
mu2 <- 2

# reguular univ. grid on 1st dimension
x1 <- seq(-2,5, len=100)
head(x1)

# reguular univ. grid on 2nd dimension
x2 <- seq(-1, 6, len=100)
head(x2)

# regular bivariate grid expanded row-wise
x1x2 <- expand.grid(x1, x2)
head(x1x2)


## ---------------------------------------------------------------------------------------
sigma1.2 <- 1
sigma2.2 <- 2
sigma12 <- 1
Sigma <- matrix(c(1,1,1,2),ncol=2)


## ---------------------------------------------------------------------------------------
z.pdf <- apply(x1x2,MARGIN = 1,
               function(x) dbvnorm(x,
                                   mu1=mu1,
                                   mu2=mu2,
                                   sigma1.2=sigma1.2,
                                   sigma2.2=sigma2.2,
                                   sigma12=sigma12))
# we reorder the object to  100 x 100 matrix
z.pdf.mat <- matrix(z.pdf, ncol=100, nrow=100, byrow = F)


## ---------------------------------------------------------------------------------------
persp(x1,x2,z.pdf.mat, zlab = "f(x)",
      col="lightblue",theta=45,phi=30,
      xlab = "x1", ylab = "x2", shade=0.01, border = 2)
contour(x1,x2,z.pdf.mat, xlab="x1", ylab="x2")
