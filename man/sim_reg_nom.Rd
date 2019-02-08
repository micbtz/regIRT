\name{sim_reg_nom}
\alias{sim_reg_nom}
\title{Simulations}
\description{
Simulation of dataset and estimation of the parameters with penalized likelihood.
Tuning parameter selected by cross-validation.
}
\usage{
sim_reg_nom(r, n, param0, lambda, itmsel, cat0, nodes, weights, pen, adaptive)
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
\item{pen}{"lasso" or "ridge"}
\item{adaptive}{logical}
}
