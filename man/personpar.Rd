\name{personpar}
\alias{personpar}
\alias{personpar.raschmodel}
\alias{personpar.rsmodel}
\alias{personpar.pcmodel}
\alias{personpar.nplmodel}
\alias{personpar.gpcmodel}

\alias{coef.personpar}
\alias{print.personpar}
\alias{vcov.personpar}

\title{Extract Person Parameters of Item Response Models}

\description{
  A class and generic function for representing and estimating the person
  parameters of a given item response model.
}

\usage{
personpar(object, \dots)
\method{personpar}{raschmodel}(object, personwise = FALSE, ref = NULL,
  vcov = TRUE, interval = NULL, tol = 1e-8, \dots)
\method{personpar}{rsmodel}(object, personwise = FALSE, ref = NULL,
  vcov = TRUE, interval = NULL, tol = 1e-8, \dots)
\method{personpar}{pcmodel}(object, personwise = FALSE, ref = NULL,
  vcov = TRUE, interval = NULL, tol = 1e-8, \dots)
\method{personpar}{nplmodel}(object, personwise = FALSE, vcov = TRUE,
  interval = NULL, tol = 1e-6, method = "EAP", \dots)
\method{personpar}{gpcmodel}(object, personwise = FALSE, vcov = TRUE,
  interval = NULL, tol = 1e-6, method = "EAP", \dots)
}

\arguments{
  \item{object}{a fitted model object for which person parameters should be
    returned/estimated.}
  \item{personwise}{logical. Should the distributional parameters of the
    latent person ability distribution be computed (default) or the person-wise
    (individual) person parameters? See below for details.}
  \item{ref}{a vector of labels or position indices of item parameters or a
    contrast matrix which should be used as restriction/for normalization. This
    argument will be passed over to internal calls of \code{itempar}.}
  \item{vcov}{logical. Should a covariance matrix be returned/estimated for the
    person parameter estimates? See also details below.}
  \item{interval}{numeric vector of length two, specifying an interval for
    \code{\link[stats]{uniroot}} or \code{\link[mirt]{fscores}} to calculate
    the person parameter estimates.}
  \item{tol}{numeric tolerance passed to \code{\link[stats]{uniroot}} or
    \code{\link[mirt]{fscores}}.}
  \item{method}{type of estimation method being passed to
    \code{\link[mirt]{fscores}}.}
  \item{\dots}{further arguments which are passed to \code{\link[stats]{optim}}
    in case of \code{vcov = TRUE} or \code{\link[mirt]{fscores}}.}
}

\details{
  \code{personpar} is both a class to represent person parameters of item
  response models as well as a generic function. The generic function can be
  used to return/estimate the person parameters of a given item response model.
  
  By default, the function \code{personpar()} reports the distribution
  parameters of the assumed person ability distribution. For models estimated by
  marginal maximum likelihood estimation (MML) this is the mean/variance of the
  underlying normal distribution, whereas for models estimated by conditional
  maximum likelihood estimation (CML) this is a discrete distribution with one
  estimation for each observed raw score in the data.
  
  Alternatively, when setting \code{personwise = TRUE}, the person parameter for
  each person/subject in the underlying data set can be extracted. In the CML
  case, this simply computes the raw score for each person and then extracts
  the corresponding person parameter. In the MML case, this necessitates
  (numerically) integrating out the individual person parameters (also known as
  factor scores or latent trait estimates) based on the underlying normal
  distribution.
  
  More specifically, the following algorithms are employed for obtaining the
  distributional person parameters: \itemize{
  
    \item In the MML case -- i.e., for \code{nplmodel}s and \code{gpcmodel}s --
      the distributional parameters are already part of the model specification.
      In a single-group specification and in the reference group of a multi-group
      specification the mean/variance parameters are fixed to 0/1. In the multi-group
      case the remaining mean/variance parameters were already estimated along with
      all other model parameters and simply need to be extracted. Analogously,
      the corresponding variance-covariance matrix just needs to be extracted and
      has zero covariances in the cells corresponding to fixed parameters.
      
    \item In the CML case -- i.e., \code{raschmodel}s, \code{rsmodel}s, and \code{pcmodel}s --
      the distributional parameters are estimated via \code{uniroot()} with the estimation
      equations given by Hoijtink & Boomsma (1995) as well as Andersen (1995). This
      approach is fast and estimates for all possible raw scores are available. If
      the covariance matrix of the estimated person parameters is requested
      (\code{vcov = TRUE}), an additional call of \code{optim} is employed to
      obtain the Hessian numerically. With this approach, person parameters are
      available only for observed raw scores.
  
  }
  
  As already explained above, obtaining the person-wise (individual) person
  paremeters (or ability estimates or factor scores) is straightforward in the
  CML case. In the MML case, \code{\link[mirt]{fscores}} is used, see Chalmers
  (2012) for further details. If \code{personwise = TRUE}, the associated
  variance-covariance matrix is not provided and simply a matrix with \code{NA}s
  is returned. (The same is done for \code{vcov = FALSE}.)
  
  For objects of class \code{personpar}, several methods to standard generic
  functions exist: \code{print}, \code{coef}, \code{vcov}. \code{coef} and
  \code{vcov} can be used to extract the person parameters and covariance matrix
  without additional attributes. Based on this Wald tests or confidence
  intervals can be easily computed, e.g., via \code{confint}.
}

\value{
  A named vector with person parmeters of class \code{personpar} and
  additional attributes \code{"model"} (the model name), \code{"vcov"} (the
  covariance matrix of the estimates if \code{vcov = TRUE} or an
  \code{NA}-matrix otherwise) and \code{"type"} (the type of the parameters,
  depending on \code{personwise}).
}

\references{
  Andersen EB (1995).
    Polytomous Rasch Models and Their Estimation.
    In Fischer GH, Molenaar IW (eds.).
    \emph{Rasch Models: Foundations, Recent Developments, and Applications}.
    Springer, New York.

  Chalmers RP (2012).
    mirt: A Multidimensional Item Response Theory Package for the R
    Environment.
    \emph{Journal of Statistical Software}, \bold{48}(6), 1--29.

  Hoijtink H, Boomsma A (1995).
    On Person Parameter Estimation in the Dichotomous Rasch Model.
    In Fischer GH, Molenaar IW (eds.).
    \emph{Rasch Models: Foundations, Recent Developments, and Applications}.
    Springer, New York.
}

\seealso{\code{\link{itempar}}, \code{\link{threshpar}},
  \code{\link{discrpar}}, \code{\link{guesspar}}, \code{\link{upperpar}}}


\examples{
o <- options(digits = 3)

## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit a Rasch model to dichotomized verbal aggression data and
ram <- raschmodel(VerbalAggression$resp2)

## extract person parameters
## (= parameters of the underlying ability distribution)
rap <- personpar(ram)
rap

## extract variance-covariance matrix and standard errors
vc <- vcov(rap)
sqrt(diag(vc))

## Wald confidence intervals
confint(rap)

## now match each person to person parameter with the corresponding raw score
personpar(ram, personwise = TRUE)[1:6]

## person parameters for RSM/PCM fitted to original polytomous data
rsm <- rsmodel(VerbalAggression$resp)
pcm <- pcmodel(VerbalAggression$resp)
cbind(personpar(rsm, vcov = FALSE), personpar(pcm, vcov = FALSE))

if(requireNamespace("mirt")) {
## fit a 2PL accounting for gender impact and
twoplm <- nplmodel(VerbalAggression$resp2, impact = VerbalAggression$gender)

## extract the person parameters
## (= mean/variance parameters from the normal ability distribution)
twoplp <- personpar(twoplm)
twoplp

## extract the standard errors
sqrt(diag(vcov(twoplp)))

## Wald confidence intervals
confint(twoplp)

## now look at the individual person parameters
## (integrated out over the normal ability distribution)
personpar(twoplm, personwise = TRUE)[1:6]
}

options(digits = o$digits)
}

\keyword{classes}

