\name{pcmodel}
\alias{pcmodel}
\alias{PCModel.fit}

\alias{print.pcmodel}
\alias{summary.pcmodel}
\alias{print.summary.pcmodel}

\alias{coef.pcmodel}
\alias{bread.pcmodel}
\alias{estfun.pcmodel}
\alias{logLik.pcmodel}
\alias{vcov.pcmodel}

\title{Partial Credit Model Fitting Function}

\description{
  \code{pcmodel} is a basic fitting function for partial credit models.
}

\usage{
pcmodel(y, weights = NULL, nullcats = c("keep", "downcode", "ignore"),
  start = NULL, reltol = 1e-10, deriv = c("sum", "diff"),
  hessian = TRUE, maxit = 100L, full = TRUE, \dots)
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
  \item{nullcats}{character string, specifying how items with
    null categories (i.e., categories not observed) should be treated (see
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
  \code{pcmodel} provides a basic fitting function for partial
  credit models, intended as a building block for fitting partial
  credit trees. It estimates the partial credit model suggested
  by Masters (1982) under the cumulative threshold parameterization,
  i.e., the item-category parameters \eqn{\eta_{jk} = \sum_{\ell =
      1}^{k}\delta_{jk}} are estimated by the the function \code{pcmodel}.
  
  Null categories, i.e., categories which have not been used, can be
  problematic when estimating a partial credit model. Several strategies
  have been suggested to cope with null categories.  \code{pcmodel}
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
  models via \code{eRm::PCM}.

  \code{pcmodel} returns an object of class \code{"pcmodel"} for
  which several basic methods are available, including \code{print},
  \code{plot}, \code{summary}, \code{coef}, \code{vcov}, \code{logLik},
  \code{\link{discrpar}}, \code{\link{itempar}}, \code{estfun},
  \code{\link{threshpar}}, and \code{\link{personpar}}.
}

\value{
  \code{pcmodel} returns an S3 object of class \code{"pcmodel"}, 
  i.e., a list the following components:
  \item{coefficients}{a named vector of estimated item-category
    parameters (without the first item-category parameter which is
    constrained to 0),}
  \item{vcov}{covariance matrix of the parameters in the model,}
  \item{data}{modified data, used for model-fitting, i.e., cleaned for
    items without variance, centralized so that the first category is
    zero for all items, treated null categories as specified via
    argument \code{"nullcats"} and without observations with zero
    weight. Be careful, this is different than for objects of class
    \code{"raschmodel"} or \code{"btmodel"}, where \code{data} contains
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

\references{
  Masters GN (1992).
    A Rasch Model for Partial Credit Scoring.
    \emph{Psychometrika}, \bold{47}(2), 149--174.

  Wilson M, Masters GN (1993).
    The Partial Credit Model and Null Categories.
    \emph{Psychometrika}, \bold{58}(1), 87--99.
}

\seealso{\code{\link{gpcmodel}}, \code{\link{rsmodel}}, \code{\link{raschmodel}},
  \code{\link{plmodel}}, \code{\link{btmodel}}}

\examples{
o <- options(digits = 4)

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Partial credit model for the other-to-blame situations
pcm <- pcmodel(VerbalAggression$resp[, 1:12])
summary(pcm)

## visualizations
plot(pcm, type = "profile")
plot(pcm, type = "regions")
plot(pcm, type = "piplot")
plot(pcm, type = "curves")
plot(pcm, type = "information")

## Get data of situation 1 ('A bus fails to
## stop for me') and induce a null category in item 2.
pcd <- VerbalAggression$resp[, 1:6, drop = FALSE]
pcd[pcd[, 2] == 1, 2] <- NA

## fit pcm to these data, comparing downcoding and keeping strategy
pcm_va_keep  <- pcmodel(pcd, nullcats = "keep")
pcm_va_down  <- pcmodel(pcd, nullcats = "downcode")

plot(x = coef(pcm_va_keep), y = coef(pcm_va_down),
     xlab = "Threshold Parameters (Keeping)",
     ylab = "Threshold Parameters (Downcoding)",
     main = "Comparison of two null category strategies (I2 with null category)", 
     pch = rep(as.character(1:6), each = 2)[-3])
abline(b = 1, a = 0)

options(digits = o$digits)
}

\keyword{regression}
