\name{nominalCV}
\alias{nominalCV}
\title{Cross-validation of a nominal model with penalty}
\description{
Selection of the tuning parameter of a nominal model with penalty using cross-validation.
}
\usage{
nominalCV(object, K, trace = FALSE)
}
\arguments{
\item{object}{output of function \code{\link{nominalmod}}.}
\item{K}{Number of folds.}
\item{trace}{Trace information.}
}
\value{
A list with components:
\item{est}{list containing the parameter estimates on the training sets. Each component corresponds to a different value of lambda, and it is a matrix with K (= number of folds) rows and J (= number of parameters) columns.}
\item{lik}{list containing the negative log-likelihood (not penalized) computed on the validation set. Each component corresponds to a different value of lambda, and it is a vector with length K.}
\item{lambda}{vector of tuning parameters.}
\item{sel}{integer value indicating the selection.}
\item{lambdasel}{value of lambda selected.}
\item{par}{matrix of parameter estimates obtained with function nominalmod. Columns correspond to different values of lambda.}
\item{data}{dataset.}
\item{D}{number of dimensions.}
}
\seealso{
\code{\link{nominalmod}}, \code{\link{regIRT}}, \code{\link{regPath}}
}




