% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ry_moebius.R
\name{ry_moebius}
\alias{ry_moebius}
\title{Mobius transformation}
\usage{
ry_moebius(zs, ws)
}
\arguments{
\item{zs}{is a complex vector with length three.}

\item{ws}{is a complex vector with length three.}
}
\description{
This function returns a moebius transformation
function that maps three complex points, zs, to
three complex points, ws.
}
\examples{
zs <- runif(3) + 1i * runif(3)
ws <- runif(3) + 1i * runif(3)
mob.fx <- ry_moebius(zs,ws)
print(mob.fx(zs))
print(ws)
}
\keyword{complex}
\keyword{conformal}
\keyword{function,}
\keyword{transformation}

