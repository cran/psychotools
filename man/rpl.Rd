\name{rpl}
\alias{rpl}

\title{Simulate Data under a Parametric Logistic IRT Model}

\description{
  \code{rpl} simulates IRT data under a parametric logistic IRT model of type
  "2PL", "3PL", "3PLu", "4PL", and "Rasch/1PL".
}

\usage{
rpl(theta, a, b, g, u, return_setting = TRUE)
}

\arguments{
  \item{theta}{numeric vector of person parameters. Can also be a list, then a
    list of length \code{length(theta)} is returned, containing multiple
    simulated data matrices.}
  \item{a}{numeric vector of item discrimination parameters.}
  \item{b}{numeric vector of item difficulty parameters.}
  \item{g}{numeric vector of so-called item guessing parameters.}
  \item{u}{numeric vector of item upper asymptote parameters.}
  \item{return_setting}{logical. Should a list containing slots of "a", "b",
    "g", "u", and "theta", as well as the simulated data matrix "data" be
    returned (default) or only the simulated data matrix.}
}


\value{
  \code{rpl} returns either a list of the following components:
  \item{a}{numeric vector of item discrimination parameters used,}
  \item{b}{numeric vector of item difficulty parameters used,}
  \item{g}{numeric vector of item guessing parameters used,}
  \item{u}{numeric vector of item upper asymptote parameters used,}
  \item{theta}{numeric vector of person parameters used,}
  \item{data}{numeric matrix containing the simulated data,}
  or (if \code{return_setting = FALSE}) only the numeric matrix containing the
  simulated data.
}

\seealso{\code{\link{rrm}}, \code{\link{rgpcm}}, \code{\link{rpcm}},
  \code{\link{rrsm}}}

\examples{
set.seed(1)
## item responses under a 2PL (two-parameter logistic) model from
## 6 persons with three different person parameters
## 9 increasingly difficult items and corresponding discrimination parameters
## no guessing (= 0) and upper asymptote 1
ppar <- rep(c(-2, 0, 2), each = 2)
ipar <- seq(-2, 2, by = 0.5)
dpar <- rep(c(0.5, 1, 1.5), each = 3)
sim <- rpl(theta = ppar, a = dpar, b = ipar, g = rep(0, 9), u = rep(1, 9))

## simulated item response data along with setting parameters
sim

## print and plot corresponding item response object
iresp <- itemresp(sim$data)
iresp
plot(iresp)
}
