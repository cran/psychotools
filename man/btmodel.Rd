\name{btmodel}
\alias{btmodel}
\alias{btReg.fit}

\alias{print.btmodel}
\alias{summary.btmodel}
\alias{print.summary.btmodel}
\alias{coef.btmodel}
\alias{worth.btmodel}
\alias{deviance.btmodel}
\alias{logLik.btmodel}
\alias{vcov.btmodel}
\alias{estfun.btmodel}


\title{Bradley-Terry Model Fitting Function}

\description{
  \code{btmodel} is a basic fitting function for simple Bradley-Terry models.
}

\usage{
btmodel(y, weights = NULL, type = c("loglin", "logit"), ref = NULL,
  undecided = NULL, position = NULL, start = NULL, vcov = TRUE, estfun =
  FALSE, \dots)
}

\arguments{
  \item{y}{paircomp object with the response.}
  \item{weights}{an optional vector of weights (interpreted as case weights).}
  \item{type}{character. Should an auxiliary log-linear Poisson model or logistic
    binomial be employed for estimation? The latter is not available if
    undecided effects are estimated.}
  \item{ref}{character or numeric. Which object parameter should be the reference
    category, i.e., constrained to zero?}
  \item{undecided}{logical. Should an undecided parameter be estimated?}
  \item{position}{logical. Should a position effect be estimated?}
  \item{start}{numeric. Starting values when calling \code{\link[stats]{glm.fit}}.}
  \item{vcov}{logical. Should the estimated variance-covariance be included in
    the fitted model object?}
  \item{estfun}{logical. Should the empirical estimating functions (score/gradient
    contributions) be included in the fitted model object?}
  \item{\dots}{further arguments passed to functions.}
}

\details{
  \code{btmodel} provides a basic fitting function for Bradley-Terry models,
  intended as a building block for fitting Bradley-Terry trees and
  Bradley-Terry mixtures in the \pkg{psychotree} package, respectively. While
  \code{btmodel} is intended for individual paired-comparison data, the
  \pkg{eba} package provides functions for aggregate data.
  
  \code{btmodel} returns an object of class \code{"btmodel"} for which
  several basic methods are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, \code{estfun}
  and \code{\link{worth}}.
}

\value{
  \code{btmodel} returns an S3 object of class \code{"btmodel"},
  i.e., a list with components as follows.
  \item{y}{paircomp object with the response}
  \item{coefficients}{estimated parameters on log-scale (without the first
    parameter which is always constrained to be 0),}
  \item{vcov}{covariance matrix of the parameters in the model,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{df}{number of estimated parameters,}
  \item{weights}{the weights used (if any),}
  \item{n}{number of observations (with non-zero weights),}
  \item{type}{character for model type (see above),}
  \item{ref}{character for reference category (see above),}
  \item{undecided}{logical for estimation of undecided parameter (see above),}
  \item{position}{logical for estimation of position effect (see above),}
  \item{labels}{character labels of the objects compared,}
  \item{estfun}{empirical estimating function (also known as scores or gradient contributions).}
}

\seealso{\code{\link{pcmodel}}, \code{\link{rsmodel}},
  \code{\link{raschmodel}}, the \pkg{eba} package}

\examples{
o <- options(digits = 4)

## data
data("GermanParties2009", package = "psychotools")

## Bradley-Terry model
bt <- btmodel(GermanParties2009$preference)
summary(bt)
plot(bt)

options(digits = o$digits)
}

\keyword{regression}
