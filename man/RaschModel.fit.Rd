\name{RaschModel.fit}
\alias{RaschModel.fit}

\alias{print.RaschModel}
\alias{summary.RaschModel}
\alias{print.summary.RaschModel}

\alias{coef.RaschModel}
\alias{worth.RaschModel}
\alias{logLik.RaschModel}
\alias{vcov.RaschModel}


\title{Rasch Model Fitting Function}

\description{
  \code{RaschModel.fit} is a basic fitting function for simple Rasch models.
}

\usage{
RaschModel.fit(y, weights = NULL, start = NULL, gradtol = 1e-6,
  deriv = c("sum", "diff", "numeric"), hessian = TRUE,
  iterlim = 100, \dots)
}

\arguments{
  \item{y}{object that can be coerced (via \code{\link[base]{as.matrix}}) to a binary 0/1 matrix.}
  \item{weights}{an optional vector of weights, interpreted as case weights (integer only).}
  \item{start}{an optional vector of starting values.}
  \item{deriv}{character. Which type of derivatives should be used for computing
    gradient and Hessian matrix? Analytical with sum algorithm (\code{"sum"}),
    analytical with difference algorithm (\code{"diff"}, faster but numerically unstable),
    or numerical.}
  \item{hessian}{logical. Should the Hessian of the final model be computed?
    If set to \code{FALSE}, the \code{vcov} method can only return \code{NA}s
    and consequently no standard errors or tests are available in the \code{summary}.}
  \item{gradtol, iterlim, \dots}{further arguments passed to \code{\link[stats]{nlm}}.}
}

\details{
  \code{RaschModel.fit} provides a basic fitting function for simple Rasch models,
  intended as a building block for fitting Rasch trees and Rasch mixtures
  in the \pkg{psychotree} package, respectively.
  
  \code{RaschModel.fit} returns an object of class \code{"RaschModel"} for which
  several basic methods are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, and \code{\link{worth}}.
}

\value{
  \code{RaschModel.fit} returns an S3 object of class \code{"RaschModel"}, 
  i.e., a list with components as follows.
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
  \item{nlm_code}{convergence code from \code{nlm},}
  \item{iterations}{number of iterations used by \code{nlm},}
  \item{gradtol}{tolerance passed to \code{nlm}.}
}

\seealso{\code{\link{btReg.fit}}}

\examples{
## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Rasch model for the self-to-blame situations
m <- RaschModel.fit(VerbalAggression$resp2[, 1:12])
summary(m)
plot(m)
}

\keyword{regression}
