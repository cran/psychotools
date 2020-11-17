\name{rrsm}
\alias{rrsm}

\title{Simulate Data under a Rating Scale Model}

\description{
  \code{rrsm} simulates IRT data under a rating scale model.}

\usage{
rrsm(theta, beta, tau, nullcats = FALSE, return_setting = TRUE)
}

\arguments{
  \item{theta}{numeric vector of person parameters. Can also be a list, then a
    list of length \code{length(theta)} is returned, containing multiple
    simulated data sets.}
  \item{beta}{numeric vector of item difficulty parameters.}
  \item{tau}{numeric vector of threshold parameters.}
  \item{nullcats}{logical. Should null categories be allowed?}
  \item{return_setting}{logical. Should a list containing slots of "beta",
    "tau", and "theta", as well as the simulated data matrix "data" be
    returned (default) or only the simulated data matrix?}
}

\value{
  \code{rrsm} returns either a list of the following components:
  \item{beta}{numeric vector of item difficulty parameters used,}
  \item{tau}{numeric vector of threshold parameters used,}
  \item{theta}{numeric vector (or list) of person parameters used,}
  \item{data}{numeric matrix containing the simulated data,}
  or (if \code{return_setting = FALSE}) only the numeric matrix containing the
  simulated data.
}

\seealso{\code{\link{rpcm}}, \code{\link{rgpcm}}, \code{\link{rrm}},
  \code{\link{rpl}}}

\examples{
set.seed(1)
## item responses under a rating scale model (RSM) with
## 6 persons with three different person parameters
## 9 increasingly difficult items
## 3 different threshold parameters
ppar <- rep(-1:1, each = 2)
ipar <- seq(-2, 2, by = 0.5)
tpar <- 0:2
sim <- rrsm(theta = ppar, beta = ipar, tau = tpar)

## simulated item response data along with setting parameters
sim

## print and plot corresponding item response object
iresp <- itemresp(sim$data)
iresp
plot(iresp)
}
