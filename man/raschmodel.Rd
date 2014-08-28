\name{raschmodel}
\alias{raschmodel}
\alias{RaschModel.fit}

\alias{print.raschmodel}
\alias{summary.raschmodel}
\alias{print.summary.raschmodel}

\alias{coef.raschmodel}
\alias{bread.raschmodel}
\alias{estfun.raschmodel}
\alias{logLik.raschmodel}
\alias{vcov.raschmodel}

\title{Rasch Model Fitting Function}

\description{
  \code{raschmodel} is a basic fitting function for simple Rasch models.
}

\usage{
raschmodel(y, weights = NULL, start = NULL, reltol = 1e-10,
  deriv = c("sum", "diff", "numeric"), hessian = TRUE,
  maxit = 100L, full = TRUE, gradtol = reltol, iterlim = maxit, \dots)
}

\arguments{
  \item{y}{object that can be coerced (via \code{\link[base]{as.matrix}}) to a binary 0/1 matrix.}
  \item{weights}{an optional vector of weights (interpreted as case weights).}
  \item{start}{an optional vector of starting values.}
  \item{deriv}{character. Which type of derivatives should be used for computing
    gradient and Hessian matrix? Analytical with sum algorithm (\code{"sum"}),
    analytical with difference algorithm (\code{"diff"}, faster but numerically unstable),
    or numerical.}
  \item{hessian}{logical. Should the Hessian of the final model be computed?
    If set to \code{FALSE}, the \code{vcov} method can only return \code{NA}s
    and consequently no standard errors or tests are available in the \code{summary}.}
  \item{reltol, maxit, \dots}{further arguments passed to \code{\link[stats]{optim}}.}
  \item{full}{logical. Should a full model object be returned? If set to \code{FALSE},
    no variance-covariance matrix and no matrix of estimating functions are computed.}
  \item{gradtol, iterlim}{numeric. For backward compatibility with previous versions
    these arguments are mapped to \code{reltol} and \code{maxit}, respectively.}
}

\details{
  \code{raschmodel} provides a basic fitting function for simple Rasch models,
  intended as a building block for fitting Rasch trees and Rasch mixtures
  in the \pkg{psychotree} and \pkg{psychomix} packages, respectively.
  
  \code{raschmodel} returns an object of class \code{"raschmodel"} for which
  several basic methods are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, \code{estfun},
  \code{\link{discrpar}}, \code{\link{itempar}}, \code{\link{threshpar}},
  and \code{\link{personpar}}.
}

\value{
  \code{raschmodel} returns an S3 object of class \code{"raschmodel"}, 
  i.e., a list with the following components:
  \item{coefficients}{estimated item difficulty parameters (without first item
    parameter which is always constrained to be 0),}
  \item{vcov}{covariance matrix of the parameters in the model,}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{df}{number of estimated parameters,}
  \item{data}{the original data supplied (excluding columns without variance),}
  \item{weights}{the weights used (if any),}
  \item{n}{number of observations (with non-zero weights),}
  \item{items}{status indicator (0, 0/1, 1) of all original items,}
  \item{na}{logical indicating whether the data contains NAs,}
  \item{elementary_symmetric_functions}{List of elementary symmetric functions
    for estimated parameters (up to order 2; or 1 in case of numeric derivatives),}
  \item{code}{convergence code from \code{optim},}
  \item{iterations}{number of iterations used by \code{optim},}
  \item{reltol}{tolerance passed to \code{optim},}
  \item{deriv}{type of derivatives used for computing gradient and
  Hessian matrix.}
}

\seealso{\code{\link{pcmodel}}, \code{\link{rsmodel}}, \code{\link{btmodel}}}

\examples{
o <- options(digits = 4)

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Rasch model for the other-to-blame situations
m <- raschmodel(VerbalAggression$resp2[, 1:12])
summary(m)

## visualizations
plot(m, type = "profile")
plot(m, type = "regions")
plot(m, type = "curves")
plot(m, type = "information")
plot(m, type = "piplot")

options(digits = o$digits)
}

\keyword{regression}
