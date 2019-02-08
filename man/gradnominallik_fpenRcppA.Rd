\name{gradnominallik_fpenRcppA}
\alias{gradnominallik_fpenRcppA}
\title{gradient of the log-likelihood of a nominal model with lasso-type penalty}
\description{
gradient of the log-likelihood of a nominal model with lasso-type penalty on slope parameters
}
\usage{
gradnominallik_fpenRcppA(par, data, nodes, weights, numpatt, itemsselect, w, 
lambda = 0, eps = 0.001)
}
\arguments{
\item{par}{vector of parameters}
\item{data}{dataset (patterns not repeated)}
\item{nodes}{quadrature nodes}
\item{weights}{quadrature weights}
\item{numpatt}{frequency of patterns}
\item{itemsselect}{items with penalty}
\item{w}{weights for adaptive penalization}
\item{lambda}{tuning parameter}
\item{eps}{constant used to approximate the absolute value}
}
