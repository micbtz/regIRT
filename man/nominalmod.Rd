\name{nominalmod}
\alias{nominalmod}
\title{Estimation of a nominal response model}
\description{
Estimates the parameters of the nominal response model with optional penalty on the slope parameters.
}
\usage{
nominalmod(data, par, lambda = 0, pen = NULL, adaptive = NULL, 
  items.select = 1:ncol(data), nq = 61)
}
\arguments{
\item{data}{Dataset.}
\item{par}{Initial values for the parameters.}
\item{lambda}{Vector of tuning parameters.}
\item{pen}{Type of penalization: "lasso" or "ridge".}
\item{adaptive}{logical; if TRUE adaptive lasso is performed.}
\item{items.select}{Items with penalty.}
\item{nq}{Number of quadrature points.}
}
\details{If lambda is zero, no penalitazion is applied.
If lambda contains more elements, the model is estimated for each value.}



