\name{regIRT-package}
\alias{regIRT-package}
\alias{regIRT}
\docType{package}
\title{Regularization methods for IRT models}
\description{Estimation of the nominal response model using a penalized likelihood approach.}
\keyword{ package }
\examples{
library(regIRT)
data(nomdata)

# estimate the nominal response model
library(mirt)
mod <- mirt(nomdata, 1, 'nominal')
par <- coef(mod)

# reparameterize the parameters
param <- lapply(par, FUN = repar)
param <- param[1:(length(param) - 1)]
par <- unlist(param)

# small lasso penalty
mod_lasso1 <- nominalmod(data = nomdata, D = 1, parini = par, lambda = 10^-3, pen = "lasso", 
adaptive = FALSE)
par_lasso1 <- mod_lasso1$par

# values of the tuning parameter lambda
lambda <- seq(0.1, 1, l = 10)

# non adaptive penalty
mods_nonadp <- nominalmod(data = nomdata, D = 1, parini = par_lasso1, lambda = lambda, 
pen = "lasso", adaptive = FALSE)
mods_nonadp$par # parameter estimates
# cross-validation
cvres_nonadp <- nominalCV(mods_nonadp, K = 5)
# regularization path
par(mfrow = c(2, 2))
regPath(cvres_nonadp)

# adaptive penalty
mods_adp <- nominalmod(data = nomdata, D = 1, parini = par_lasso1, parW = par_lasso1, 
lambda = lambda, pen = "lasso", adaptive = TRUE)
# cross-validation
cvres_adp <- nominalCV(mods_adp, K = 5)
# regularization path
par(mfrow = c(2, 2))
regPath(cvres_adp)

}

