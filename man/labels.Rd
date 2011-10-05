\name{labels<-}

\alias{labels<-}

\title{Set Labels}

\description{
  A generic function for setting labels for an object.
}

\usage{
  labels(object) <- value
}

\arguments{
  \item{object}{an object.}
  \item{value}{an object.}
}

\examples{
## method for "paircomp" data
pc <- paircomp(rbind(
  c(1,  1,  1), # a > b, a > c, b > c
  c(1,  1, -1), # a > b, a > c, b < c
  c(1, -1, -1), # a > b, a < c, b < c
  c(1,  1,  1)))
labels(pc)
labels(pc) <- c("ah", "be", "ce")
pc
}

\keyword{classes}
