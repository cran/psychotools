\name{personpar}
\alias{personpar}
\alias{personpar.raschmodel}
\alias{personpar.rsmodel}
\alias{personpar.pcmodel}

\alias{coef.personpar}
\alias{print.personpar}
\alias{vcov.personpar}

\title{Extract Person Parameters of Item Response Models}

\description{
  A class and generic function for representing and estimating
  the person parameters of a given item response model.
}

\usage{
personpar(object, \dots)
\method{personpar}{raschmodel}(object, ref = NULL, vcov = TRUE,
  interval = NULL, tol = 1e-8, \dots)
\method{personpar}{rsmodel}(object, ref = NULL, vcov = TRUE,
  interval = NULL, tol = 1e-8, \dots)
\method{personpar}{pcmodel}(object, ref = NULL, vcov = TRUE,
  interval = NULL, tol = 1e-8, \dots)
}

\arguments{
  \item{object}{a fitted model object for which person parameters
    should be estimated.}
  \item{ref}{a vector of labels or position indices of item parameters
    or a contrast matrix which should be used as restriction/for
    normalization. If \code{NULL} (the default), all items are used
    (sum zero restriction). This argument will passed over to
    internal calls of \code{itempar}.}
  \item{vcov}{logical. Should a covariance matrix be estimated (based
    on the joint log-likelihood) for the person parameter estimates?
    See also details below.}
  \item{interval}{numeric vector of length two, specifying an interval
    for \code{\link[stats]{uniroot}} to search for the person parameter estimates.}
  \item{tol}{numeric tolerance passed to \code{\link[stats]{uniroot}}.}
  \item{\dots}{further arguments which are passed to \code{\link[stats]{optim}}
    in case of \code{vcov = TRUE}.}
}

\details{
  \code{personpar} is both, a class to represent person parameters
  of item response models as well as a generic function. The generic
  function can be used to estimate the person parameters of a given item
  response model.

  Person parameters are estimated via \code{uniroot()} with
  the estimation equations given by Hoijtink & Boomsma (1995) as well
  as Andersen (1995). This approach is fast and estimates for all
  possible raw scores are available. If the covariance matrix of the
  estimated person parameters is requested (\code{vcov = TRUE}), an
  additional call of \code{optim} is necessary to obtain the Hessian.
  With this approach, person parameters are available only for observed
  raw scores.

  For objects of class \code{personpar}, several methods to standard
  generic functions exist: \code{print}, \code{coef}, \code{vcov}.
  \code{coef} and \code{vcov} can be used to extract the
  person parameters and covariance matrix without additional
  attributes. Based on this Wald tests or confidence intervals can be
  easily computed, e.g., via \code{confint}.
}

\value{
  A named vector with estimated person parmeters of class
  \code{personpar} and additional attributes \code{"model"} (the model
  name), and \code{"vcov"} (the covariance matrix of the estimates if
  \code{vcov = TRUE} or an NA-matrix otherwise).
}

\references{
  Erling B. Andersen. Polytomous Rasch Models and Their Estimation. In
  Gerhard H. Fischer, and Ivo W. Molenaar, (1995). Rasch Models:
  Foundations, Recent Developments, and Applications.

  Herbert Hoijtink, and Anne Boomsma. On Person Parameter Estimation in
  the Dichotomous Rasch Model. In Gerhard H. Fischer, and Ivo
  W. Molenaar, (1995). Rasch Models: Foundations, Recent Developments,
  and Applications.  
}

\seealso{\code{\link{itempar}}, \code{\link{threshpar}},
  \code{\link{discrpar}}}

\examples{
o <- options(digits = 4)

## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit a Rasch model to dichotomized verbal aggression data and
ram <- raschmodel(VerbalAggression$resp2)

## extract person parameters
rap <- personpar(ram)
rap

## extract variance-covariance matrix and standard errors
vc <- vcov(rap)
sqrt(diag(vc))

## Wald confidence intervals
confint(rap)

## person parameters for RSM/PCM fitted to original polytomous data
rsm <- rsmodel(VerbalAggression$resp)
pcm <- pcmodel(VerbalAggression$resp)
cbind(personpar(rsm, vcov = FALSE), personpar(pcm, vcov = FALSE))

options(digits = o$digits)
}

\keyword{classes}

