\name{covariates}

\alias{covariates}
\alias{covariates<-}

\title{Extract/Set Covariates}

\description{
  A generic function for extracting/setting covariates for an object.
}

\usage{
  covariates(object, \dots)
  covariates(object) <- value
}

\arguments{
  \item{object}{an object.}
  \item{\dots}{arguments passed to methods.}
  \item{value}{an object.}
}

\examples{
## method for "paircomp" data
pc <- paircomp(rbind(
  c(1,  1,  1), # a > b, a > c, b > c
  c(1,  1, -1), # a > b, a > c, b < c
  c(1, -1, -1), # a > b, a < c, b < c
  c(1,  1,  1)))
covariates(pc)
covariates(pc) <- data.frame(foo = factor(c(1, 2, 2), labels = c("foo", "bar")))
covariates(pc)
}

\keyword{classes}
