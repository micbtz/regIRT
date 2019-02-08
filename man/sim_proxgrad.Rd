\name{sim_proxgrad}
\alias{sim_proxgrad}
\title{Simulations}
\description{
Simulation of dataset and estimation of the parameters with penalized likelihood
using the proximal gradient algorithm.
}
\usage{
sim_proxgrad(r, n, param0, lambda, itmsel, cat0, nodes, weights, adaptive)
}
\arguments{
\item{r}{replication number}
\item{n}{sample size}
\item{param0}{true item parameters}
\item{lambda}{tuning parameter}
\item{itmsel}{items with penalty}
\item{cat0}{number of categories for each item}
\item{nodes}{quadrature nodes}
\item{weights}{quadrature weights}
\item{adaptive}{logical}
}
