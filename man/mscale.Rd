\name{mscale}

\alias{mscale}
\alias{mscale<-}

\title{Extract/Replace Measurement Scale}

\description{
  Generic functions for extracting and replacing the measurement scale from an object.
}

\usage{
  mscale(object, \dots)
  mscale(object) <- value
}

\arguments{
  \item{object}{an object.}
  \item{\dots}{arguments passed to methods.}
  \item{value}{an object describing the measurement scale.}
}

\examples{
## methods for "paircomp" data
pc <- paircomp(rbind(
  c(2,  1,  0),
  c(1,  1, -1),
  c(1, -2, -1),
  c(0,  0,  0)))
pc

## extract
mscale(pc)

## replace (collapse to >/=/< scale)
mscale(pc) <- sign(mscale(pc))
pc
}

\keyword{classes}
