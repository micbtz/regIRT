\name{nominalmod}
\alias{nominalmod}
\title{Estimation of a nominal response model}
\description{
Estimates the parameters of the nominal response model with optional penalty on the slope parameters.
}
\usage{
nominalmod(data, D, parini, parW = NULL, lambda = 0, pen = NULL, adaptive = NULL, 
  items.select = 1:ncol(data), nq = NULL)
}
\arguments{
\item{data}{dataset.}
\item{D}{number of dimensions.}
\item{parini}{initial values for the parameters.}
\item{parW}{vector of parameters used for computing weights for the adaptive penalization.}
\item{lambda}{vector of tuning parameters.}
\item{pen}{type of penalization: "lasso" or "ridge".}
\item{adaptive}{logical; if TRUE adaptive lasso is performed.}
\item{items.select}{vector of integer values indicating the items with penalty.}
\item{nq}{number of quadrature points per dimension. By default the number of quadrature points depends on the number of dimensions: '1'=61, '2'=31, '3'=15, '4'=9, '5'=7, '>5'=3.}
}
\details{If lambda is zero, no penalitazion is applied.
If lambda contains more elements, the model is estimated for each value.}
\value{
A list with components:
\item{data}{dataset.}
\item{D}{number of dimensions.}
\item{parini}{initial values for the parameters.}
\item{parW}{vector of parameters used for computing weights for the adaptive penalization.}
\item{lambda}{vector of tuning parameters.}
\item{pen}{type of penalization: "lasso" or "ridge".}
\item{adaptive}{logical; if TRUE adaptive lasso is performed.}
\item{items.select}{vector of integer values indicating the items with penalty.}
\item{nq}{number of quadrature points per dimension.}
\item{par}{matrix of parameter estimates. Columns correspond to different values of lambda.}
\item{lik}{vector containing the penalized log-likelihood computed for each lambda.}
\item{convergence}{An integer code. 0 indicates successful completion.}
}
\seealso{
\code{\link{nominalCV}}, \code{\link{regIRT}}, \code{\link{regPath}}
}




