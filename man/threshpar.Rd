\name{threshpar}
\alias{threshpar}
\alias{threshpar.raschmodel}
\alias{threshpar.rsmodel}
\alias{threshpar.pcmodel}

\alias{coef.threshpar}
\alias{print.threshpar}

\title{Extract Threshold Parameters of Item Response Models}

\description{
  A class and generic function for representing and extracting
  the item threshold parameters of a given item response model.
}

\usage{
  threshpar(object, \dots)
  \method{threshpar}{raschmodel}(object, type = c("mode", "median",
  "mean"), ref = NULL, alias = TRUE, relative = FALSE,
   cumulative = FALSE, vcov = TRUE, \dots)
  \method{threshpar}{rsmodel}(object, type = c("mode", "median",
  "mean"), ref = NULL, alias = TRUE, relative = FALSE,
   cumulative = FALSE, vcov = TRUE, \dots)
  \method{threshpar}{pcmodel}(object, type = c("mode", "median",
   "mean"), ref = NULL, alias = TRUE, relative = FALSE,
   cumulative = FALSE, vcov = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object whose threshold parameters
    should be extracted.}
  \item{type}{character of length one which determines the type of
  threshold parameters to return (see details below).}
  \item{ref}{a vector of labels or position indices of (relative)
    threshold parameters or a contrast matrix which should be used
    as restriction/for normalization. For partial credit models,
    argument \code{ref} can also be a list of m contrasts. If
    \code{NULL} (the default), the relative threshold parameters
    are centered around their item-specific means and the absolute
    threshold parameters are centered around their global mean.}
  \item{alias}{logical. If \code{TRUE} (the default), the aliased
    parameter is included in the return vector (and in the variance-covariance
    matrix if \code{vcov} = TRUE). If \code{FALSE}, it is removed. If the
    restriction given in \code{ref} depends on several parameters, the
    first parameter of the restriction specified is (arbitrarily) chosen
    to be removed if \code{alias} is \code{FALSE}.}
  \item{relative}{logical. If set to \code{FALSE} (default), absolute item
    threshold parameters are returned. If set to \code{TRUE}, relative
    item threshold parameters with the contrast specified
    in argument \code{ref} are returned.}
  \item{cumulative}{logical. If set to \code{TRUE}, cumulative threshold
    parameters are returned. These correspond to the cumulative sum over
    the absolute or relative item threshold parameters (after the
    restriction given in argument \code{ref} has been applied).}
  \item{vcov}{logical. If \code{TRUE} (the default), the (transformed)
    variance-covariance matrix of the (relative) threshold parameters
    is attached as attribute \code{vcov}.. If \code{FALSE}, a NA-matrix
    is attached.}
  \item{\dots}{further arguments which are currently not used.}
}

\details{
  \code{threshpar} is both, a class to represent threshold parameters
  of item response models as well as a generic function. The generic
  function can be used to extract the threshold parameters of a given
  item response model.

  For objects of class \code{threshpar}, methods to standard generic
  functions \code{print} and \code{coef} can be used to print and
  extract the threshold parameters.

  Depending on argument \code{type}, different item threshold parameters
  are returned. For \code{type = "mode"}, the returned item threshold
  parameters correspond to the location on the theta axis where the
  probability of category \eqn{k} equals the probability of category
  \eqn{k-1}. For Rasch and partial credit models, item threshold parameters
  of this type correspond directly to the estimated absolute item
  threshold parameters of these models. For \code{type = "median"}, the
  returned item threshold parameters correspond to the location on the
  theta axis where the probability of choosing category \eqn{k} or higher,
  i.e., \eqn{P(X_{ij} >= k)}, equals 0.5. For \code{type = "mean"}, the
  returned absolute item threshold parameters correspond to the location
  on the theta axis where the expected category response is in the
  middle between two categories, i.e. 0.5, 1.5, \dots{}. An illustration
  of these threshold parameters can be found on page 104 in Masters &
  Wright (1995).
}

\value{
  A named list with item threshold parameters of class \code{threshpar} and
  additional attributes \code{model} (the model name), \code{type}
  (the type of item threshold parameters returned, see details above),
  \code{ref} (the items or parameters used as restriction/for
  normalization), \code{relative} (whether relative or absolute
  item threshold parameters are returned), \code{cumulative} (whether
  the cumulative item threshold parameters are returned),
  \code{alias} (either \code{FALSE} or a named character vector or list
  with the removed aliased parameters), and \code{vcov} (the estimated
  and adjusted variance-covariance matrix).
}

\references{
  Masters GN, Wright BD (1997).
    The Partial Credit Model.
    In Van der Linden WJ, Hambleton RK (eds.).
    \emph{Handbook of Modern Item Response Theory}.
}

\seealso{\code{\link{itempar}}, \code{\link{personpar}}, and \code{\link{discrpar}}}

\examples{
o <- options(digits = 4)

## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit a rasch model to dichotomized verbal aggression data
raschmod <- raschmodel(VerbalAggression$resp2)

## extract threshold parameters with sum zero restriction
tr <- threshpar(raschmod)
tr

## compare to item parameters (again with sum zero restriction)
ip <- itempar(raschmod)
ip

all.equal(coef(tr), coef(ip))

## rating scale model example
rsmod <- rsmodel(VerbalAggression$resp)
trmod <- threshpar(rsmod, type = "mode")
trmed <- threshpar(rsmod, type = "median")
trmn <- threshpar(rsmod, type = "mean")

## compare different types of threshold parameters
cbind("Mode" = coef(trmod, type = "vector"),
 "Median" = coef(trmod, type = "vector"),
 "Mean" = coef(trmn, type = "vector"))

options(digits = o$digits)
}

\keyword{classes}
