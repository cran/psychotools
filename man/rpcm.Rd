\name{rpcm}
\alias{rpcm}

\title{Simulate Data under a Partial Credit Model}

\description{
  \code{rpcm} simulates IRT data under a partial credit model.}

\usage{
rpcm(theta, delta, nullcats = FALSE, return_setting = TRUE)
}

\arguments{
  \item{theta}{numeric vector of person parameters. Can also be a list, then a
    list of length \code{length(theta)} is returned, containing multiple
    simulated data sets.}
  \item{delta}{list of numeric vectors of item threshold parameters.}
  \item{nullcats}{logical. Should null categories be allowed?}
  \item{return_setting}{logical. Should a list containing slots of "delta",
    and "theta", as well as the simulated data matrix "data" be
    returned (default) or only the simulated data matrix?}
}

\value{
  \code{rpcm} returns either a list of the following components:
  \item{delta}{list of numeric vectors of item threshold parameters used,}
  \item{theta}{numeric vector of person parameters used,}
  \item{data}{numeric matrix containing the simulated data,}
  or (if \code{return_setting = FALSE}) only the numeric matrix containing the
  simulated data.
}

\seealso{\code{\link{rgpcm}}, \code{\link{rrsm}}, \code{\link{rrm}},
  \code{\link{rpl}}}

\examples{
set.seed(1)
## item responses under a partial credit model (PCM) with
## 6 persons with three different person parameters
## 8 items with different combinations of two or three threshold parameters
ppar <- rep(-1:1, each = 2)
tpar <- rep(list(-2:0, -1:1, 0:1, 0:2), each = 2)
sim <- rpcm(theta = ppar, delta = tpar)

## simulated item response data along with setting parameters
sim

## print and plot corresponding item response object
iresp <- itemresp(sim$data)
iresp
plot(iresp)
}
