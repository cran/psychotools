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

\seealso{\code{\link{btReg.fit}}, \code{\link{RaschModel.fit}}}

\examples{
## data
data("GermanParties2009", package = "psychotools")

## Bradley-Terry model
bt <- btReg.fit(GermanParties2009$preference)

## worth parameters
worth(bt)
}

\keyword{classes}
