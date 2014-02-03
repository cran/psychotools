\name{RSModel.fit}
\alias{RSModel.fit}

\alias{print.RSModel}
\alias{summary.RSModel}
\alias{print.summary.RSModel}

\alias{coef.RSModel}
\alias{logLik.RSModel}
\alias{vcov.RSModel}


\title{Rating Scale Model Fitting Function}

\description{

  \code{RSModel.fit} is a basic fitting function for rating scale models.

}


\usage{

RSModel.fit(y, weights = NULL, start = NULL, reltol = 1e-10,
  deriv = c("sum", "diff"), hessian = TRUE,
  maxit = 100L, full = TRUE, \dots)

}

\arguments{

  \item{y}{object that can be coerced (via
  \code{\link[base]{as.matrix}}). Typically either already a
  \code{\link[base]{matrix}} or a \code{\link[base]{data.frame}} with
  items in the columns and observations in the rows.}

  \item{weights}{an optional vector of weights, interpreted as case
    weights (integer only).}
  
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

  \code{RSModel.fit} provides a basic fitting function for rating scales models,
  intended as a building block for fitting rating scale trees.
  
  \code{RSModel.fit} returns an object of class \code{"RSModel"} (and
  class \code{"PCModel"}) for which several basic methods
  are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, and \code{logLik}.

}

\value{

  \code{RSModel.fit} returns an S3 object of class \code{"RSModel"}, 
  i.e., a list with components as follows.

  \item{coefficients}{a named vector of estimated item difficulty
    parameters (without first item parameter which is constrained
    to 0) and estimated treshold parameters (again without first
    threshold parameter which is also constrained to 0),}
  
  \item{vcov}{covariance matrix of the parameters in the model,}
  
  \item{data}{modified data, used for model-fitting, i.e., cleaned for items without
    variance, centralized so that the first category is zero for all items
    and without observations with zero weight. Be careful, this is different than for
    objects of class \code{"RaschModel"} or \code{"btReg"}, where
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

  \item{reltol}{tolerance passed to \code{optim}.}

}


\seealso{\code{\link{PCModel.fit}}, \code{\link{RaschModel.fit}}, \code{\link{btReg.fit}}}


\examples{

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Rating scale model for the other-to-blame situations
rsm <- RSModel.fit(VerbalAggression$resp[, 1:12])
summary(rsm)
plot(rsm)

}


\keyword{regression}
