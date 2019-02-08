\name{nominalCV}
\alias{nominalCV}
\title{Cross-validation of a nominal model with penalty}
\description{
Cross-validation of a nominal model with penalty
}
\usage{
nominalCV(data, K, par, lambda, nodes, weights, items.select = 1:ncol(data), 
pen, trace = FALSE, w)

}
\arguments{
\item{data}{dataset}
\item{K}{number of folds}
\item{par}{vector of parameters}
\item{lambda}{tuning parameter}
\item{nodes}{quadrature nodes}
\item{weights}{quadrature weights}
\item{items.select}{items with penalty}
\item{pen}{"lasso" or "ridge"}
\item{trace}{trace information}
\item{w}{weights for adaptive penalty}
}
