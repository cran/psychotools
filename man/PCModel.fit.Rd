\name{PCModel.fit}
\alias{PCModel.fit}

\alias{print.PCModel}
\alias{summary.PCModel}
\alias{print.summary.PCModel}

\alias{coef.PCModel}
\alias{logLik.PCModel}
\alias{vcov.PCModel}


\title{Partial Credit Model Fitting Function}

\description{
  \code{PCModel.fit} is a basic fitting function for partial credit models.
}

\usage{
PCModel.fit(y, weights = NULL, nullcats = c("keep", "downcode", "ignore"),
  start = NULL, reltol = 1e-10, deriv = c("sum", "diff"),
  hessian = TRUE, maxit = 100L, full = TRUE, \dots)
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

  \item{nullcats}{character vector of length one, specifying how items with
  null categories, i.e., categories not used, should be treated (see
  details below).}

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

  \code{PCModel.fit} provides a basic fitting function for partial credit models,
  intended as a building block for fitting partial credit trees.
  
  Null categories, i.e., categories which have not been used, can be
  problematic when estimating a partial credit model. Several strategies
  have been suggested to cope with null categories.  \code{PCModel.fit}
  allows to select from three possible strategies via the argument
  \code{nullcats}. If \code{nullcats} is set to \code{"keep"} (the
  default), the strategy suggested by Wilson & Masters (1993) is used to
  handle null categories. That basically means that the integrity of the
  response framework is maintained, i.e., no category scores are
  changed. This is not the case, when \code{nullcats} is set to
  \code{"downcode"}. Then all categories above a null category are
  shifted down to close the existing gap. In both cases (\code{"keep"}
  and \code{"downcode"}) the number of estimated parameters is reduced
  by the number of null categories. When \code{nullcats} is set to
  \code{"ignore"}, these are literally ignored and a threshold parameter
  is estimated during the optimization nevertheless. This strategy is
  used by the related package \pkg{eRm} when fitting partial credit
  models via \code{\link[eRm]{PCM}}.

  \code{PCModel.fit} returns an object of class \code{"PCModel"} for
  which several basic methods are available, including \code{print},
  \code{plot}, \code{summary}, \code{coef}, \code{vcov}, and \code{logLik}.

}

\value{

  \code{PCModel.fit} returns an S3 object of class \code{"PCModel"}, 
  i.e., a list with components as follows.

  \item{coefficients}{a nemd vector of estimated item-category
    parameters (without first item-category parameter which is
    constrained to 0),}
  
  \item{vcov}{covariance matrix of the parameters in the model,}
  
  \item{data}{modified data, used for model-fitting, i.e., cleaned for
    items without variance, centralized so that the first category is
    zero for all items, treated null categories as specified via
    argument \code{"nullcats"} and without observations with zero
    weight. Be careful, this is different than for objects of class
    \code{"RaschModel"} or \code{"btReg"}, where \code{data} contains
    the \emph{original} data,}

  \item{items}{logical vector of length \code{ncol(dat)}, indicating
    which items have variance (\code{TRUE}), i.e., are identified and
    have been used for the estimation or not (\code{FALSE}),}
    
  \item{categories}{list of length \code{ncol(y)}, containing integer
    vectors starting from one to the number of categories minus one per
    item,}

  \item{n}{number of observations (with non-zero weights),}

  \item{n_org}{original number of observations in \code{y},}

  \item{weights}{the weights used (if any),}
  
  \item{na}{logical indicating whether the data contain NAs,}

  \item{nullcats}{either \code{NULL} or, if there have been null
  categories, a list of length \code{ncol(y)} with logical vectors
  specifying which categories are null categories (\code{TRUE}) or not
  (\code{FALSE}),}

  \item{esf}{list of elementary symmetric functions and their
    derivatives for estimated parameters,}

  \item{loglik}{log-likelihood of the fitted model,}

  \item{df}{number of estimated parameters,}

  \item{code}{convergence code from \code{optim},}

  \item{iterations}{number of iterations used by \code{optim},}

  \item{reltol}{tolerance passed to \code{optim}.}

}

\references{Wilson, M. & Masters, G. N. (1993). The partial credit model
  and null categories. Psychometrika, 58(1), 87-99.}

\seealso{\code{\link{RSModel.fit}}, \code{\link{RaschModel.fit}}, \code{\link{btReg.fit}}}

\examples{

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Partial credit model for the other-to-blame situations
pcm <- PCModel.fit(VerbalAggression$resp[, 1:12])
summary(pcm)
plot(pcm)

## Get data of situation 1 ('A bus fails to
## stop for me') and induce a null category in item 2.
pcd <- VerbalAggression$resp[, 1:6, drop = FALSE]
pcd[pcd[, 2] == 1, 2] <- NA

## fit pcm to these data, comparing downcoding and keeping strategy
pcm_va_keep  <- PCModel.fit(pcd, nullcats = "keep")
pcm_va_down  <- PCModel.fit(pcd, nullcats = "downcode")

plot(x = coef(pcm_va_keep), y = coef(pcm_va_down),
     xlab = "Threshold Parameters (Downcoding)",
     ylab = "Threshold Parameters (Keeping)",
     main = "Comparison of two null category strategies (I2 with null category)", 
     pch = rep(as.character(1:6), each = 2)[-3])
abline(b = 1, a = 0)

}

\keyword{regression}
