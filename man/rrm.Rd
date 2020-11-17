\name{rrm}
\alias{rrm}

\title{Simulate Data under a Rasch model}

\description{
  \code{rrm} simulates IRT data under a Rasch model.
}

\usage{
rrm(theta, beta, return_setting = TRUE)
}

\arguments{
  \item{theta}{numeric vector of person parameters. Can also be a list, then a
    list of length \code{length(theta)} is returned, containing multiple
    simulated data matrices.}
  \item{beta}{numeric vector of item difficulty parameters.}
  \item{return_setting}{logical. Should a list containing slots of "beta",
    and "theta", as well as the simulated data matrix "data" be
    returned (default) or only the simulated data matrix.}
}


\value{
  \code{rrm} returns either a list of the following components:
  \item{beta}{numeric vector of item difficulty parameters used,}
  \item{theta}{numeric vector of person parameters used,}
  \item{data}{numeric matrix containing the simulated data,}
  or (if \code{return_setting = FALSE}) only the numeric matrix containing the
  simulated data.
}

\seealso{\code{\link{rpl}}, \code{\link{rpcm}}, \code{\link{rrsm}},
  \code{\link{rgpcm}}}

\examples{
set.seed(1)
## item responses under a Rasch model from
## 6 persons with three different person parameters
## 9 increasingly difficult items
ppar <- rep(-1:1, each = 2)
ipar <- seq(-2, 2, by = 0.5)
sim <- rrm(theta = ppar, beta = ipar)

## simulated item response data along with setting parameters
sim

## print and plot corresponding item response object
iresp <- itemresp(sim$data)
iresp
plot(iresp)
}
