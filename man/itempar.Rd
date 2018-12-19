\name{itempar}
\alias{itempar}
\alias{itempar.btmodel}
\alias{itempar.raschmodel}
\alias{itempar.rsmodel}
\alias{itempar.pcmodel}
\alias{itempar.plmodel}
\alias{itempar.gpcmodel}
\alias{itempar.raschtree}
\alias{itempar.bttree}

\alias{coef.itempar}
\alias{print.itempar}
\alias{vcov.itempar}

\title{Extract Item Parameters of Item Response Models}

\description{
  A class and generic function for representing and extracting the item
  parameters of a given item response model.
}

\usage{
  itempar(object, \dots)
  \method{itempar}{raschmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{itempar}{rsmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{itempar}{pcmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{itempar}{plmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{itempar}{gpcmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{itempar}{btmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, log = FALSE, \dots)
}

\arguments{
  \item{object}{a fitted model or tree object whose item parameters should be extracted.}
  \item{ref}{a vector of labels or position indices of item parameters or a
    contrast matrix which should be used as restriction/for normalization. If
    \code{NULL} (the default) for all models except models estimated via MML,
    all items are used (sum zero restriction). For models estimated via MML
    (\code{plmodel}s and \code{gpcmodel}s), the parameters are by default
    identified via the distributional parameters of the person parameters (mean
    and variance of a normal distribution). Nevertheless, a restriction on the
    interval scale can be applied.}
  \item{alias}{logical. If \code{TRUE} (the default), the aliased parameter is
    included in the return vector (and in the variance-covariance matrix if
    \code{vcov} = TRUE). If \code{FALSE}, it is removed. If the restriction
    given in \code{ref} depends on several parameters, the first parameter of
    the restriction specified is (arbitrarily) chosen to be removed if
    \code{alias} is \code{FALSE}.}
  \item{vcov}{logical. If \code{TRUE} (the default), the (transformed)
    variance-covariance matrix of the item parameters is attached as
    attribute \code{vcov}. If \code{FALSE}, an \code{NA}-matrix is attached.}
  \item{log}{logical. Whether to return the estimated model parameters
    on the logit (\code{TRUE}) or preference scale (\code{FALSE}).}
  \item{\dots}{further arguments which are currently not used.}
}

\details{
  \code{itempar} is both, a class to represent item parameters of item
  response models as well as a generic function. The generic function can be
  used to extract the item parameters of a given item response model.

  For Rasch models and n-parameter logistic models, \code{itempar} returns the
  estimated item difficulty parameters \eqn{\hat{\beta}_{j}} under the
  restriction specified in argument \code{ref}. For rating scale models,
  \code{itempar} returns computed item location parameters \eqn{\hat{\beta}_{j}}
  under the restriction specified in argument \code{ref}. These are computed
  from the estimated item-specific parameters \eqn{\hat{\xi}_{j}} (who mark the
  location of the first category of an item on the latent theta axis). For
  partial credit models and generalized partial credit models, \code{itempar}
  returns \sQuote{mean} absolute item threshold parameters, \eqn{\hat{\beta}_{j}
  = \frac{1}{p_{j}} \sum_{k = 1}^{p_{j}}\hat{\delta}_{jk}}, i.e., a single
  parameter per item is returned which results as the mean of the absolute item
  threshold parameters \eqn{\hat{\delta}_{jk}} of this item. Based upon these
  \sQuote{mean} absolute item threshold parameters \eqn{\hat{\beta}_{j}}, the
  restriction specified in argument \code{ref} is applied. For all models, the
  variance-covariance matrix of the returned item parameters is adjusted
  according to the multivariate delta rule.

  For objects of class \code{itempar}, several methods to standard generic
  functions exist: \code{print}, \code{coef}, \code{vcov}. \code{coef} and
  \code{vcov} can be used to extract the estimated calculated item parameters
  and their variance-covariance matrix without additional attributes. Based on
  this Wald tests or confidence intervals can be easily computed, e.g., via
  \code{confint}.

  Two-sample item-wise Wald tests for DIF in the item parameters can be
  carried out using the function \code{\link{anchortest}}.
}

\value{
  A named vector with item parameters of class \code{itempar} and additional
  attributes \code{model} (the model name), \code{ref} (the items or parameters
  used as restriction/for normalization), \code{alias} (either \code{FALSE} or a
  named character vector with the removed aliased parameter, and \code{vcov}
  (the adjusted covariance matrix of the estimates if \code{vcov = TRUE} or an
  \code{NA}-matrix otherwise).
}

\seealso{\code{\link{personpar}}, \code{\link{threshpar}},
  \code{\link{discrpar}}, \code{\link{guesspar}}, \code{\link{upperpar}}}

\examples{
o <- options(digits = 4)

## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit a Rasch model to dichotomized verbal aggression data
raschmod <- raschmodel(VerbalAggression$resp2)

## extract item parameters with sum zero or use last two items as anchor
ip1 <- itempar(raschmod)
ip2a <- itempar(raschmod, ref = 23:24) # with position indices
ip2b <- itempar(raschmod, ref = c("S4WantShout", "S4DoShout")) # with item label

ip1
ip2a

all.equal(ip2a, ip2b)

## extract vcov
vc1 <- vcov(ip1)
vc2 <- vcov(ip2a)

## adjusted standard errors,
## smaller with more items used as anchors
sqrt(diag(vc1))
sqrt(diag(vc2))

## Wald confidence intervals
confint(ip1)
confint(ip2a)

options(digits = o$digits)
}

\keyword{classes}

