\name{nominal_grfused}
\alias{nominal_grfused}
\title{ADMM algorithm}
\description{
Estimation of the parameters with penalized likelihood
using the ADMM algorithm.
}
\usage{
nominal_grfused(par, data, maxiter = 1000, maxiter_innerloop = 100, v = NULL, 
cinit, lambda, items.select = 1:ncol(data), nodes, weights)
}
\arguments{
\item{par}{initial values}
\item{data}{dataset}
\item{maxiter}{maximum number of iterations}
\item{maxiter_innerloop}{maximum number of iterations of the inner loop}
\item{v}{intial values for lagrange multipliers}
\item{cinit}{initial value for the constant that multiplies the quadratic term}
\item{lambda}{tuning parameter}
\item{items.select}{items with penalty}
\item{nodes}{quadrature nodes}
\item{weights}{quadrature weights}
}
