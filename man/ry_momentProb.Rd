% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ry_momentProb.R
\name{ry_momentProb}
\alias{ry_momentProb}
\title{Arbitrary moment}
\usage{
ry_momentProb(x, p, order = 1, center = FALSE)
}
\arguments{
\item{x}{is a numeric vector.}

\item{p}{is a numeric vector of weight of x elements}

\item{order}{is order of moment.}

\item{center}{is logical; when TRUE, moment is around mean and otherwise zero.}
}
\description{
This function returns order-th weighted sample moment.
}
\examples{
n <- 10
x <- runif(n)
p <- runif(n)
p <- p/sum(p)
out <- ry_momentProb(x,p,oder=2,center=TRUE)
}
\keyword{moment}

