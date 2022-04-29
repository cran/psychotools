\name{discrpar}
\alias{discrpar}
\alias{discrpar.raschmodel}
\alias{discrpar.rsmodel}
\alias{discrpar.pcmodel}
\alias{discrpar.nplmodel}
\alias{discrpar.gpcmodel}

\alias{coef.discrpar}
\alias{print.discrpar}
\alias{vcov.discrpar}

\title{Extract Discrimination Parameters of Item Response Models}

\description{
  A class and generic function for representing and extracting the
  discrimination parameters of a given item response model.
}

\usage{
  discrpar(object, \dots)
  \method{discrpar}{raschmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{discrpar}{rsmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{discrpar}{pcmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{discrpar}{nplmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{discrpar}{gpcmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object whose discrimination parameters should be
    extracted.}
  \item{ref}{a restriction to be used. Not used for models estimated via CML as
    the discrimination parameters are fixed to 1 in \code{raschmodel}s,
    \code{rsmodel}s and \code{pcmodel}s. For models estimated via MML
    (\code{nplmodel}s and \code{gpcmodel}s), the parameters are by default
    identified via the distributional parameters of the person parameters (mean and
    variance of a normal distribution). Nevertheless, a restriction on the ratio
    scale can be applied.}
  \item{alias}{logical. If \code{TRUE} (the default), the aliased parameters are
    included in the return vector (and in the variance-covariance matrix if
    \code{vcov} = TRUE). If \code{FALSE}, these parameters are removed. For
    \code{raschmodel}s, \code{rsmodel}s and \code{pcmodel}s where all
    discrimination parameters are fixed to 1, this means that an empty
    numeric vector and an empty variance-covariance matrix is returned if
    \code{alias} is \code{FALSE}.}
  \item{vcov}{logical. If \code{TRUE} (the default), the variance-covariance
    matrix of the discrimination parameters is attached as attribute
    \code{vcov}.}
  \item{\dots}{further arguments which are currently not used.}
}

\details{
  \code{discrpar} is both, a class to represent discrimination parameters of
  item response models as well as a generic function. The generic function can
  be used to extract the discrimination parameters of a given item response
  model.

  For objects of class \code{discrpar}, several methods to standard generic
  functions exist: \code{print}, \code{coef}, \code{vcov}. \code{coef} and
  \code{vcov} can be used to extract the discrimination parameters and their
  variance-covariance matrix without additional attributes.
}

\value{
  A named vector with discrimination parameters of class \code{discrpar} and
  additional attributes \code{model} (the model name), \code{ref} (the items or
  parameters used as restriction/for normalization), \code{alias} (either
  \code{TRUE} or a named numeric vector with the aliased parameters not included
  in the return value), and \code{vcov} (the estimated and adjusted
  variance-covariance matrix).
}

\seealso{\code{\link{personpar}}, \code{\link{itempar}},
  \code{\link{threshpar}}, \code{\link{guesspar}}, \code{\link{upperpar}}}

\examples{
o <- options(digits = 4)

## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit Rasch model to verbal aggression data
rmod <- raschmodel(VerbalAggression$resp2)

## extract the discrimination parameters
dp1 <- discrpar(rmod)

## extract the standard errors
sqrt(diag(vcov(dp1)))

if(requireNamespace("mirt")) {
## fit 2PL to verbal aggression data
twoplmod <- nplmodel(VerbalAggression$resp2)

## extract the discrimination parameters
dp2 <- discrpar(twoplmod)

## this time with the first discrimination parameter being the reference
discrpar(twoplmod, ref = 1)

## extract the standard errors
sqrt(diag(vcov(dp2)))
}

options(digits = o$digits)
}

\keyword{classes}

