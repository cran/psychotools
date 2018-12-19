\name{upperpar}
\alias{upperpar}
\alias{upperpar.raschmodel}
\alias{upperpar.rsmodel}
\alias{upperpar.pcmodel}
\alias{upperpar.plmodel}
\alias{upperpar.gpcmodel}

\alias{coef.upperpar}
\alias{print.upperpar}
\alias{vcov.upperpar}

\title{Extract Upper Asymptote Parameters of Item Response Models}

\description{
  A class and generic function for representing and extracting the
  upper asymptote parameters of a given item response model.
}

\usage{
  upperpar(object, \dots)
  \method{upperpar}{raschmodel}(object, alias = TRUE, vcov = TRUE, \dots)
  \method{upperpar}{rsmodel}(object, alias = TRUE, vcov = TRUE, \dots)
  \method{upperpar}{pcmodel}(object, alias = TRUE, vcov = TRUE, \dots)
  \method{upperpar}{plmodel}(object, alias = TRUE, logit = FALSE, vcov = TRUE, \dots)
  \method{upperpar}{gpcmodel}(object, alias = TRUE, vcov = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object whose upper asymptote parameters should be
    extracted.}
  \item{alias}{logical. If \code{TRUE} (the default), the aliased parameters
    are included in the return vector (and in the variance-covariance matrix
    if \code{vcov} = TRUE). If \code{FALSE}, these parameters are removed. For
    \code{raschmodel}s, \code{rsmodel}s, \code{pcmodel}s and \code{gpcmodel}s,
    where all upper asymptote parameters are fixed to 1, this means
    that an empty numeric vector and an empty variance-covariace matrix is
    returned if \code{alias} is \code{FALSE}.}
  \item{logit}{logical. If a \code{plmodel} of \code{type} \code{"3PLu"} or
    \code{"4PL"} model has been fit, the upper asymptote parameters were estimated
    on the logit scale. If \code{logit = FALSE}, these estimates and the
    variance-covariance (if requested) are retransformed using the logistic
    function and the delta method.}
  \item{vcov}{logical. If \code{TRUE} (the default), the variance-covariance
    matrix of the upper asymptote parameters is attached as attribute
    \code{vcov}.}
  \item{\dots}{further arguments which are currently not used.}
}

\details{
  \code{upperpar} is both, a class to represent upper asymptote parameters of
  item response models as well as a generic function. The generic function can
  be used to extract the upper asymptote parameters of a given item response
  model.

  For objects of class \code{upperpar}, several methods to standard generic
  functions exist: \code{print}, \code{coef}, \code{vcov}. \code{coef} and
  \code{vcov} can be used to extract the upper asymptote parameters and their
  variance-covariance matrix without additional attributes.
}

\value{
  A named vector with upper asymptote parameters of class \code{upperpar} and
  additional attributes \code{model} (the model name), \code{alias} (either
  \code{TRUE} or a named numeric vector with the aliased parameters not included
  in the return value), \code{logit} (indicating whether the estimates are on the
  logit scale or not), and \code{vcov} (the estimated and adjusted
  variance-covariance matrix).
}

\seealso{\code{\link{personpar}}, \code{\link{itempar}},
  \code{\link{threshpar}}, \code{\link{discrpar}}, \code{\link{guesspar}}}

\examples{
if(requireNamespace("mirt")) {

o <- options(digits = 3)

## load simulated data
data("Sim3PL", package = "psychotools")

## fit 2PL to data simulated under the 3PLu
twoplmod <- plmodel(Sim3PL$resp2)

## extract the upper asymptote parameters (all fixed at 1)
up1 <- upperpar(twoplmod)

## fit 3PLu to data simulated under the 3PLu
threeplmodu <- plmodel(Sim3PL$resp2, type = "3PLu")

## extract the upper asymptote parameters
up2 <- upperpar(threeplmodu)

## extract the standard errors
sqrt(diag(vcov(up2)))

## extract the upper asymptote parameters on the logit scale
up2_logit <- upperpar(threeplmodu, logit = TRUE)

## along with the delta transformed standard errors
sqrt(diag(vcov(up2_logit)))

options(digits = o$digits)
}
}

\keyword{classes}

