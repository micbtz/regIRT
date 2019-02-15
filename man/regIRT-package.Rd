\name{regIRT-package}
\alias{regIRT-package}
\alias{regIRT}
\docType{package}
\title{Regularization methods for IRT models}
\description{Estimation of the nominal response model using a penalized likelihood approach.}
\author{
	Michela Battauz
	
	Maintainer: Michela Battauz <michela.battauz@uniud.it>
}
\keyword{ package }
\examples{
data(nomdata)
# estimate the nominal response model
library(mirt)
mod <- mirt(nomdata, 1, 'nominal')
par <- coef(mod)
# reparameterize the parameters
param <- lapply(par, FUN = repar)
param <- param[1:4]
par <- unlist(param)
# collapse patterns
rdata <- reduce_data(nomdata)
datared <- rdata$data
numpatt <- rdata$numpatt
# weights for adaptive lasso
nitems <- ncol(nomdata)
w <- list()
for (j in 1:nitems) {
  m <- length(param[[j]]) / 2
  alpha <- c(0, param[[j]][1:m])
  w[[j]] <- 1 / abs(outer(alpha, alpha, "-"))
}
# quadrature points and nodes
library(statmod)
nq <- 61 # number of quadrature points
gq <- gauss.quad.prob(nq, dist = "normal")
nodes <- gq$nodes
weights <- gq$weights
# penalized estimation
lambda <- 2 #tuning parameter
opt <- optim(par = par, fn = nominallik_fpenRcppA, gr = gradnominallik_fpenRcppA,
  data = datared, method = "BFGS", lambda = lambda, nodes = nodes, weights = weights,
  numpatt = numpatt, itemsselect = 0:3, control = list(maxit = 500), eps = 0.01, w = w)
for (i in 3:9) opt <- optim(par = opt$par, fn = nominallik_fpenRcppA, 
  gr = gradnominallik_fpenRcppA, data = datared, method = "BFGS", 
  lambda = lambda, nodes = nodes, weights = weights, numpatt = numpatt, 
  itemsselect = 0:3, eps = 10^-i, w = w)
opt
}

