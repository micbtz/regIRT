\name{nominalCV}
\alias{nominalCV}
\title{Cross-validation of a nominal model with penalty}
\description{
Cross-validation of a nominal model with penalty
}
\usage{
nominalCV(data, K, par, lambda, pen, adaptive, items.select = 1:ncol(data), nq,
  trace = FALSE)
}
\arguments{
\item{data}{Dataset.}
\item{K}{Number of folds.}
\item{par}{Vector of parameters.}
\item{lambda}{Tuning parameter.}
\item{pen}{Type of penalization: "lasso" or "ridge".}
\item{adaptive}{Logical; if TRUE adaptive lasso is performed.}
\item{items.select}{Items with penalty.}
\item{nq}{Number of quadrature points.}
\item{trace}{Trace information.}
}
