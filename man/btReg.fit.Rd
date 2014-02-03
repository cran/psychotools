\name{btReg.fit}
\alias{btReg.fit}

\alias{print.btReg}
\alias{summary.btReg}
\alias{print.summary.btReg}
\alias{coef.btReg}
\alias{worth.btReg}
\alias{deviance.btReg}
\alias{logLik.btReg}
\alias{vcov.btReg}


\title{Bradley-Terry Model Fitting Function}

\description{
  \code{btReg.fit} is a basic fitting function for simple Bradley-Terry models.
}

\usage{
btReg.fit(y, weights = NULL, type = c("loglin", "logit"), ref = NULL,
  undecided = NULL, position = NULL, start = NULL, vcov = TRUE, estfun = TRUE,
  \dots)
}

\arguments{
  \item{y}{paircomp object with the response.}
  \item{weights}{an optional vector of weights, interpreted as case weights (integer only).}
  \item{type}{character. Should an auxiliary log-linear Poisson model or logistic
    binomial be employed for estimation? The latter is only available if not undecided
    effects are estimated.}
  \item{ref}{character or numeric. Which object parameter should be the reference
    category, i.e., constrained to zero?}
  \item{undecided}{logical. Should an undecided parameter be estimated?}
  \item{position}{logical. Should a position effect be estimated?}
  \item{start}{numeric. Starting values when calling \code{\link[stats]{glm.fit}}.}
  \item{vcov}{logical. Should the estimated variance-covariance be included in
    the fitted model object?}
  \item{estfun}{logical. Should the empirical estimating functions (score/gradient
    contributions be included in the fitted model object?}
  \item{\dots}{further arguments passed to functions.}
}

\details{
  \code{btReg.fit} provides a basic fitting function for Bradley-Terry models,
  intended as a building block for fitting Bradley-Terry trees and Bradley-Terry mixtures
  in the \pkg{psychotree} package, respectively.
  
  \code{btReg.fit} returns an object of class \code{"btReg"} for which
  several basic methods are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, and \code{\link{worth}}.
}

\value{
  \code{btReg.fit} returns an S3 object of class \code{"btReg"},
  i.e., a list with components as follows.
  \item{coefficients}{estimated parameters on log-scale (without the first
    parameter which is always constrained to be 0),}
  \item{vcov}{covariance matrix of the parameters in the model,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{df}{number of estimated parameters,}
  \item{estfun}{empirical estimating function (also known as scores or gradient contributions),}
  \item{weights}{the weights used (if any),}
  \item{n}{number of observations (with non-zero weights),}
  \item{type}{character for model type (see above),}
  \item{ref}{character for reference category (see above),}
  \item{undecided}{logical for estimation of undecided parameter (see above),}
  \item{position}{logical for estimation of position effect (see above),}
  \item{labels}{character labels of the objects compared.}
}

\seealso{\code{\link{PCModel.fit}}, \code{\link{RSModel.fit}}, \code{\link{RaschModel.fit}}}

\examples{
## data
data("GermanParties2009", package = "psychotools")

## Bradley-Terry model
bt <- btReg.fit(GermanParties2009$preference)
summary(bt)
plot(bt)
}

\keyword{regression}
