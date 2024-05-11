\name{rsmodel}
\alias{rsmodel}
\alias{RSModel.fit}

\alias{print.rsmodel}
\alias{summary.rsmodel}
\alias{print.summary.rsmodel}

\alias{coef.rsmodel}
\alias{bread.rsmodel}
\alias{estfun.rsmodel}
\alias{logLik.rsmodel}
\alias{vcov.rsmodel}

\title{Rating Scale Model Fitting Function}

\description{
  \code{rsmodel} is a basic fitting function for rating scale models.
}

\usage{
rsmodel(y, weights = NULL, start = NULL, reltol = 1e-10,
  deriv = c("sum", "diff"), hessian = TRUE,
  maxit = 100L, full = TRUE, \dots)
}

\arguments{
  \item{y}{item response object that can be coerced (via \code{\link[base]{as.matrix}})
    to a numeric matrix with scores 0, 1, \dots Typically, either
    already a matrix, data frame, or dedicated object of class
    \code{\link{itemresp}}.}
  \item{weights}{an optional vector of weights (interpreted as case
    weights).}
  \item{deriv}{character. If "sum" (the default), the first derivatives
    of the elementary symmetric functions are calculated with the sum
    algorithm. Otherwise ("diff") the difference algorithm (faster but
    numerically unstable) is used.}
  \item{start}{an optional vector of starting values.}
  \item{hessian}{logical. Should the Hessian of the final model be computed?
    If set to \code{FALSE}, the \code{vcov} method can only return \code{NA}s
    and consequently no standard errors or tests are available in the
    \code{summary}.}
  \item{reltol, maxit, \dots}{further arguments passed to \code{\link[stats]{optim}}.}
  \item{full}{logical. Should a full model object be returned? If set to \code{FALSE},
    no variance-covariance matrix and no matrix of estimating functions are computed.}
}

\details{
  \code{rsmodel} provides a basic fitting function for rating scales models,
  intended as a building block for fitting rating scale trees. It
  estimates the rating scale model in the parametrization suggested by
  Andrich (1978), i.e., item-specific parameters \eqn{\xi_{j}} who mark
  the location of the first absolute threshold of an item on the theta axis and
  cumulative relative threshold parameters \eqn{\kappa_{k}} are
  estimated by the function \code{rsmodel}.

  \code{rsmodel} returns an object of class \code{"rsmodel"} (and
  class \code{"pcmodel"}) for which several basic methods are available,
  including \code{print}, \code{plot}, \code{summary}, \code{coef},
  \code{vcov}, \code{logLik}, \code{\link{discrpar}}, \code{estfun},
  \code{\link{itempar}}, \code{\link{threshpar}}, and \code{\link{personpar}}.
}

\value{
  \code{rsmodel} returns an S3 object of class \code{"rsmodel"}, 
  i.e., a list with the following components:
  \item{coefficients}{a named vector of estimated item-specific
    parameters (without the first item parameter which is constrained
    to 0) and estimated cumulative relative treshold parameters
    (again without first threshold parameter which is also constrained to 0),}
  \item{vcov}{covariance matrix of the parameters in the model,}
  \item{data}{modified data, used for model-fitting, i.e., cleaned for items without
    variance, centralized so that the first category is zero for all items
    and without observations with zero weight. Be careful, this is different than for
    objects of class \code{"raschmodel"} or \code{"btmodel"}, where
    \code{data} contains the \emph{original} data,}
  \item{items}{logical vector of length \code{ncol(y)}, which
    indicates which items have variance (\code{TRUE}), i.e., are identified and have been
    used for the estimation or not (\code{FALSE}),}
  \item{categories}{integer vector of length \code{ncol(y)}, which
    contains the number of categories minus one per item,}
  \item{n}{number of observations (with non-zero weights),}
  \item{n_org}{original number of observations in \code{y},}
  \item{weights}{the weights used (if any),}
  \item{na}{logical indicating whether the data contains NAs,}
  \item{esf}{list of elementary symmetric functions and their
    derivatives for estimated parameters,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{df}{number of estimated parameters,}
  \item{code}{convergence code from \code{optim},}
  \item{iterations}{number of iterations used by \code{optim},}
  \item{reltol}{tolerance passed to \code{optim},}
  \item{call}{original function call.}
}

\references{
  Andrich D (1978).
    Application of a Psychometric Rating Model to Ordered Categories Which Are Scored with Successive Integers.  
    \emph{Psychometrika}, \bold{2}(4), 581--594.
}

\seealso{\code{\link{pcmodel}}, \code{\link{gpcmodel}}, \code{\link{raschmodel}},
  \code{\link{nplmodel}}, \code{\link{btmodel}}}

\examples{
o <- options(digits = 4)

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Rating scale model for the other-to-blame situations
rsm <- rsmodel(VerbalAggression$resp[, 1:12])
summary(rsm)

## visualizations
plot(rsm, type = "profile")
plot(rsm, type = "regions")
plot(rsm, type = "curves")
plot(rsm, type = "information")
plot(rsm, type = "piplot")

options(digits = o$digits)
}

\keyword{regression}
