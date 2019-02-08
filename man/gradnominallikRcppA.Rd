\name{gradnominallikRcppA}
\alias{gradnominallikRcppA}
\title{gradient of the log-likelihood of a nominal model}
\description{
gradient of the log-likelihood of a nominal model
}
\usage{
gradnominallikRcppA(par, data, nodes, weights, numpatt, itemsselect, lambda = 0)
}
\arguments{
\item{par}{vector of parameters}
\item{data}{dataset (patterns not repeated)}
\item{nodes}{quadrature nodes}
\item{weights}{quadrature weights}
\item{numpatt}{frequency of patterns}
\item{itemsselect}{items with penalty}
\item{lambda}{tuning parameter}
}
\details{
If lambda>0 a ridge-type penalty is added to the log-likelihood
}