\name{mptmodel}
\alias{mptmodel}

\alias{coef.mptmodel}
\alias{confint.mptmodel}
\alias{deviance.mptmodel}
\alias{estfun.mptmodel}
\alias{logLik.mptmodel}
\alias{predict.mptmodel}
\alias{print.mptmodel}
\alias{summary.mptmodel}
\alias{print.summary.mptmodel}
\alias{vcov.mptmodel}

\alias{mptspec}
\alias{print.mptspec}
\alias{update.mptspec}


\title{Multinomial Processing Tree (MPT) Model Fitting Function}

\description{
  \code{mptmodel} is a basic fitting function for multinomial processing tree
    (MPT) models.
}

\usage{
mptmodel(y, weights = NULL, spec, treeid = NULL,
  optimargs = list(control = list(reltol = .Machine$double.eps^(1/1.2),
                                  maxit = 1000),
                   init = NULL),
  start = NULL, vcov = TRUE, estfun = FALSE, \dots)
}

\arguments{
  \item{y}{matrix of response frequencies.}
  \item{weights}{an optional vector of weights (interpreted as case weights).}
  \item{spec}{an object of class \code{mptspec}: typically result of a call to
    \code{\link{mptspec}}. A symbolic description of the model to be fitted.}
  \item{treeid}{a vector that identifies each tree in a joint multinomial
    model.}
  \item{optimargs}{a list of arguments passed to the optimization function
    (\code{\link{optim}}).}
  \item{start}{a vector of starting values for the parameter estimates between
    zero and one.}
  \item{vcov}{logical. Should the estimated variance-covariance be included in
    the fitted model object?}
  \item{estfun}{logical. Should the empirical estimating functions
    (score/gradient contributions) be included in the fitted model object?}
  \item{\dots}{further arguments passed to functions.}
}

\details{
  \code{mptmodel} provides a basic fitting function for multinomial processing
  tree (MPT) models, intended as a building block for fitting MPT trees in the
  \pkg{psychotree} package. While \code{mptmodel} is intended for individual
  response frequencies, the \pkg{mpt} package provides functions for aggregate
  data.

  MPT models are specified using the \code{mptspec} function. See the
  documentation in the \pkg{mpt} package for details.
  
  \code{mptmodel} returns an object of class \code{"mptmodel"} for which
  several basic methods are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, \code{estfun}
  and \code{\link{predict}}.
}

\value{
  \code{mptmodel} returns an S3 object of class \code{"mptmodel"},
    i.e., a list with components as follows:
  \item{y}{a matrix with the response frequencies,}
  \item{coefficients}{estimated parameters (for extraction, the \code{coef}
    function is preferred),}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{npar}{number of estimated parameters,}
  \item{weights}{the weights used (if any),}
  \item{nobs}{number of observations (with non-zero weights),}
  \item{ysum}{the aggregate response frequencies,}
  \item{fitted, goodness.of.fit, ...}{see \code{mpt} in the \pkg{mpt}
    package.}
}

\seealso{\code{\link{btmodel}}, \code{\link{pcmodel}}, \code{\link{rsmodel}},
  \code{\link{raschmodel}}, \code{\link{mptspec}}, the \pkg{mpt} package}

\examples{
o <- options(digits = 4)

## data
data("SourceMonitoring", package = "psychotools")

## source-monitoring MPT model
mpt1 <- mptmodel(SourceMonitoring$y, spec = mptspec("SourceMon"))
summary(mpt1)
plot(mpt1)

options(digits = o$digits)
}

\keyword{regression}
