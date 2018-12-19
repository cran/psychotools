\name{guesspar}
\alias{guesspar}
\alias{guesspar.raschmodel}
\alias{guesspar.rsmodel}
\alias{guesspar.pcmodel}
\alias{guesspar.plmodel}
\alias{guesspar.gpcmodel}

\alias{coef.guesspar}
\alias{print.guesspar}
\alias{vcov.guesspar}

\title{Extract Guessing Parameters of Item Response Models}

\description{
  A class and generic function for representing and extracting the
  so-called guessing parameters of a given item response model.
}

\usage{
  guesspar(object, \dots)
  \method{guesspar}{raschmodel}(object, alias = TRUE, vcov = TRUE, \dots)
  \method{guesspar}{rsmodel}(object, alias = TRUE, vcov = TRUE, \dots)
  \method{guesspar}{pcmodel}(object, alias = TRUE, vcov = TRUE, \dots)
  \method{guesspar}{plmodel}(object, alias = TRUE, logit = FALSE, vcov = TRUE, \dots)
  \method{guesspar}{gpcmodel}(object, alias = TRUE, vcov = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object whose guessing parameters should be
    extracted.}
  \item{alias}{logical. If \code{TRUE} (the default), the aliased parameters
    are included in the return vector (and in the variance-covariance matrix if
    \code{vcov} = TRUE). If \code{FALSE}, these parameters are removed. For
    \code{raschmodel}s, \code{rsmodel}s, \code{pcmodel}s and \code{gpcmodel}s,
    where all guessing parameters are fixed to 0, this means that an
    empty numeric vector and an empty variance-covariace matrix is returned if
    \code{alias} is \code{FALSE}.}
  \item{logit}{logical. If a \code{plmodel} of \code{type} \code{"3PL"} or
    \code{"4PL"} model has been fit, the guessing parameters were estimated on the
    logit scale. If \code{logit = FALSE}, these estimates and the
    variance-covariance (if requested) are retransformed using the logistic
    function and the delta method.}
  \item{vcov}{logical. If \code{TRUE} (the default), the variance-covariance
    matrix of the guessing parameters is attached as attribute
    \code{vcov}.}
  \item{\dots}{further arguments which are currently not used.}
}

\details{
  \code{guesspar} is both, a class to represent guessing parameters of item
  response models as well as a generic function. The generic function can be
  used to extract the guessing parameters of a given item response model.

  For objects of class \code{guesspar}, several methods to standard generic
  functions exist: \code{print}, \code{coef}, \code{vcov}. \code{coef} and
  \code{vcov} can be used to extract the guessing parameters and their
  variance-covariance matrix without additional attributes.
}

\value{
  A named vector with guessing parameters of class \code{guesspar} and
  additional attributes \code{model} (the model name), \code{alias} (either
  \code{TRUE} or a named numeric vector with the aliased parameters not included
  in the return value), \code{logit} (indicating whether the estimates are on the
  logit scale or not), and \code{vcov} (the estimated and adjusted
  variance-covariance matrix).
}

\seealso{\code{\link{personpar}}, \code{\link{itempar}},
  \code{\link{threshpar}}, \code{\link{discrpar}}, \code{\link{upperpar}}}

\examples{
if(requireNamespace("mirt")) {

o <- options(digits = 3)

## load simulated data
data("Sim3PL", package = "psychotools")

## fit 2PL to data simulated under the 3PL
twoplmod <- plmodel(Sim3PL$resp)

## extract the guessing parameters (all fixed at 0)
gp1 <- guesspar(twoplmod)

## fit 3PL to data simulated under the 3PL
threeplmod <- plmodel(Sim3PL$resp, type = "3PL")

## extract the guessing parameters
gp2 <- guesspar(threeplmod)

## extract the standard errors
sqrt(diag(vcov(gp2)))

## extract the guessing parameters on the logit scale
gp2_logit <- guesspar(threeplmod, logit = TRUE)

## along with the delta transformed standard errors
sqrt(diag(vcov(gp2_logit)))

options(digits = o$digits)
}
}

\keyword{classes}

