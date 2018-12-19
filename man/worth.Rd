\name{worth}

\alias{worth}

\title{Extract Worth Parameters}

\description{
  Generic functions for extracting worth parameters from paired comparison models.
}

\usage{
  worth(object, \dots)
}

\arguments{
  \item{object}{an object.}
  \item{\dots}{arguments passed to methods.}
}

\details{
  Since version 0.3-0, calls to \code{worth} are
  internally passed over to \code{itempar}.
}

\seealso{\code{\link{btmodel}}, \code{\link{raschmodel}}}

\examples{
o <- options(digits = 4)

## data
data("GermanParties2009", package = "psychotools")

## Bradley-Terry model
bt <- btmodel(GermanParties2009$preference)

## worth parameters
worth(bt)

## or
itempar(bt)

options(digits = o$digits)
}

\keyword{classes}
