\name{discrpar}
\alias{discrpar}
\alias{discrpar.raschmodel}
\alias{discrpar.rsmodel}
\alias{discrpar.pcmodel}

\alias{coef.discrpar}
\alias{print.discrpar}
\alias{vcov.discrpar}

\title{Extract Discrimination Parameters of Item Response Models}

\description{
  A class and generic function for representing and extracting
  the discrimination parameters of a given item response model.
}

\usage{
  discrpar(object, \dots)
  \method{discrpar}{raschmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
  \method{discrpar}{rsmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE,  \dots)
  \method{discrpar}{pcmodel}(object, ref = NULL, alias = TRUE, vcov = TRUE, \dots)
}

\arguments{
  \item{object}{a fitted model object whose discrimination parameters
    should be extracted.}
  \item{ref}{a restriction to be used to ensure parameter
    identifiability of the item discrimination parameters. Currently not
    used as the item discrimination parameters are fixed to unity in
    \code{raschmodel}'s, \code{rsmodel}'s and \code{pcmodel}'s.}
  \item{alias}{logical. If \code{TRUE} (the default), the aliased
    parameters are included in the return vector (and in the
    variance-covariance matrix if \code{vcov} = TRUE). If \code{FALSE}, these
    parameters are removed. For \code{raschmodel}'s, \code{rsmodel}'s
    and \code{pcmodel}'s where all item discrimination parameters are
    fixed to unity, this means that an empty numeric vector and an empty
    variance-covariace matrix is returned if \code{alias} is \code{FALSE}.}
  \item{vcov}{logical. If \code{TRUE} (the default), the
    variance-covariance matrix of the discrimination parameters is
    attached as attribute \code{vcov}.}
  \item{\dots}{further arguments which are currently not used.}
}

\details{
  \code{discrpar} is both, a class to represent discrimination parameters
  of item response models as well as a generic function. The generic
  function can be used to extract the discrimination parameters of a
  given item response model.

  For objects of class \code{discrpar}, several methods to standard
  generic functions exist: \code{print}, \code{coef}, \code{vcov}.
  \code{coef} and \code{vcov} can be used to extract the discrimination
  parameters and their variance-covariance matrix without additional
  attributes.}

\value{
  A named vector with discrimination parameters of class
  \code{discrpar} and additional attributes \code{model} (the model
  name), \code{ref} (the items or parameters used as restriction/for
  normalization), \code{alias} (either \code{TRUE} or a named numeric
  vector with the aliased parameters not included in the return value),
  and \code{vcov} (the estimated and adjusted variance-covariance matrix).
}

\seealso{\code{\link{personpar}}, \code{\link{itempar}},
  \code{\link{threshpar}}}

\examples{
## load verbal aggression data
data("VerbalAggression", package = "psychotools")

## fit a rasch model to dichotomized verbal aggression data
rmod <- raschmodel(VerbalAggression$resp2)

## extract the discrimination parameters
dp <- discrpar(rmod)

## extract the variance-covariance matrix
head(vcov(dp))
}

\keyword{classes}
