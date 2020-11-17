\name{rgpcm}
\alias{rgpcm}

\title{Simulate Data under a Generalized Partial Credit Model}

\description{
  \code{rgpcm} simulates IRT data under a generalized partial credit model.}

\usage{
rgpcm(theta, a, b, nullcats = FALSE, return_setting = TRUE)
}

\arguments{
  \item{theta}{numeric vector of person parameters. Can also be a list, then a
    list of length \code{length(theta)} is returned, containing multiple
    simulated data sets.}
  \item{a}{list of numerics of item discrimination parameters.}
  \item{b}{list of numeric vectors of item threshold parameters.}
  \item{nullcats}{logical. Should null categories be allowed?}
  \item{return_setting}{logical. Should a list containing slots of "a", "b",
    and "theta", as well as the simulated data matrix "data" be
    returned (default) or only the simulated data matrix?}
}

\value{
  \code{rgpcm} returns either a list of the following components:
  \item{a}{list of numerics of item discrimination parameters used,}
  \item{b}{list of numeric vectors of item threshold parameters used,}
  \item{theta}{numeric vector of person parameters used,}
  \item{data}{numeric matrix containing the simulated data,}
  or (if \code{return_setting = FALSE}) only the numeric matrix containing the
  simulated data.
}

\seealso{\code{\link{rpcm}}, \code{\link{rrsm}}, \code{\link{rpl}},
  \code{\link{rrm}}}

\examples{
set.seed(1)
## item responses under a GPCM (generalized partial credit model) from
## 6 persons with three different person parameters
## 8 items with different combinations of two or three threshold parameters
## and corresponding discrimination parameters
ppar <- rep(-1:1, each = 2)
tpar <- rep(list(-2:0, -1:1, 0:1, 0:2), each = 2)
dpar <- rep(list(1, 2), each = 4)
sim <- rgpcm(theta = ppar, a = dpar, b = tpar)

## simulated item response data along with setting parameters
sim

## print and plot corresponding item response object
iresp <- itemresp(sim$data)
iresp
plot(iresp)
}
