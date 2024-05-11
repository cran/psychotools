\name{gpcmodel}
\alias{gpcmodel}

\alias{print.gpcmodel}
\alias{summary.gpcmodel}
\alias{print.summary.gpcmodel}

\alias{coef.gpcmodel}
\alias{bread.gpcmodel}
\alias{estfun.gpcmodel}
\alias{logLik.gpcmodel}
\alias{vcov.gpcmodel}

\title{Generalized Partial Credit Model Fitting Function}

\description{
  \code{gpcmodel} is a basic fitting function for generalized partial credit
  models providing a wrapper around \code{\link[mirt]{mirt}} and
  \code{\link[mirt]{multipleGroup}} relying on marginal maximum likelihood (MML)
  estimation via the standard EM algorithm.
}

\usage{
gpcmodel(y, weights = NULL, impact = NULL, type = c("GPCM", "PCM"),
  grouppars = FALSE, vcov = TRUE, nullcats = "downcode", 
  start = NULL, method = "BFGS", maxit = 500, reltol = 1e-5, ...)
}

\arguments{
  \item{y}{item response object that can be coerced (via \code{\link[base]{as.matrix}})
    to a numeric matrix with scores 0, 1, \dots Typically, either
    already a matrix, data frame, or dedicated object of class
    \code{\link{itemresp}}.}
  \item{weights}{an optional vector of weights (interpreted as case
    weights)}.
  \item{impact}{an optional \code{factor} allowing for grouping
    the subjects (rows). If specified, a multiple-group model is fitted
    to account for impact (see details below). By default, no impact is
    modelled, i.e., a single-group model is used.}
  \item{type}{character string, specifying the type of model to be
    estimated. In addition to the default GPCM (generalized partial credit model)
    it is also possible to estimate a standard PCM (partial credit model)
    by marginal maximum likelihood (MML).}
  \item{grouppars}{logical. Should the estimated distributional group parameters
    of a multiple group model be included in the model parameters?}
  \item{vcov}{logical or character specifying the type of variance-covariance
    matrix (if any) computed for the final model. The default \code{vcov = TRUE}
    corresponds to \code{vcov = "Oakes"}, see \code{\link[mirt]{mirt}} for
    further options. If set to \code{vcov = FALSE} (or \code{vcov = "none"}),
    \code{vcov()} will return a matrix of \code{NA}s only.}
  \item{nullcats}{character string, specifying how items with
    null categories (i.e., categories not observed) should be treated. Currently
    only \code{"downcode"} is available, i.e., all categories above a null category are
    shifted down to close the observed gap(s).}
  \item{start}{an optional vector or list of starting values (see examples
    below).}
  \item{method, maxit, reltol}{control parameters for the optimizer employed by
    \code{\link[mirt]{mirt}} for the EM algorithm.}
  \item{\dots}{further arguments passed to \code{\link[mirt]{mirt}} or
    \code{\link[mirt]{multipleGroup}}, respectively.}
}

\details{
  \code{gpcmodel} provides a basic fitting function for generalized partial
  credit models (GPCMs) providing a wrapper around \code{\link[mirt]{mirt}} (and
  \code{\link[mirt]{multipleGroup}}, respectively) relying on MML estimation via
  the standard EM algorithm (Bock & Aitkin, 1981). Models are estimated under the
  slope/intercept parametrization, see e.g. Chalmers (2012). The probability of
  person \eqn{i} falling into category \eqn{x_{ij}} of item \eqn{j} out of all
  categories \eqn{p_{j}} is modelled as:
  \deqn{P(X_{ij} = x_{ij}|\theta_{i},a_{j},\boldsymbol{d_{j}}) =
  \frac{\exp{(x_{ij}a_{j}\theta_{i} + d_{jx_{ij}})}}{\displaystyle\sum_{k = 0} ^
  {p_{j}}\exp{(ka_{j}\theta_{i} + d_{jk})}}}

  Note that all \eqn{d_{j0}} are fixed at 0. A reparametrization of the
  intercepts to the classical IRT parametrization, see e.g. Muraki (1992), is
  provided via \code{\link{threshpar}}.

  If an optional \code{impact} variable is supplied, a multiple-group model of
  the following form is being fitted: Item parameters are fixed to be equal
  across the whole sample. For the first group of the \code{impact} variable the
  person parameters are fixed to follow the standard normal distribution. In the
  remaining \code{impact} groups, the distributional parameters (mean and
  variance of a normal distribution) of the person parameters are
  estimated freely. See e.g. Baker & Kim (2004, Chapter 11), Debelak & Strobl
  (2019), or Schneider et al. (2022) for further details. To improve convergence of the model fitting
  algorithm, the first level of the \code{impact} variable should always correspond
  to the largest group. If this is not the case, levels are re-ordered internally.

  If \code{grouppars} is set to \code{TRUE} the freely estimated distributional
  group parameters (if any) are returned as part of the model parameters.

  Instead of the default GPCM, a standard partial credit model (PCM) can also
  be estimated via MML by setting \code{type = "PCM"}. In this case all slopes
  are restricted to be equal across all items.

  \code{gpcmodel} returns an object of class \code{"gpcmodel"} for which
  several basic methods are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, \code{estfun},
  \code{\link{discrpar}}, \code{\link{itempar}}, \code{\link{threshpar}}, and
  \code{\link{personpar}}.
}

\value{
  \code{gpcmodel} returns an S3 object of class \code{"gpcmodel"}, 
  i.e., a list of the following components:
  \item{coefficients}{estimated model parameters in slope/intercept parametrization,}
  \item{vcov}{covariance matrix of the model parameters,}
  \item{data}{modified data, used for model-fitting, i.e., centralized so
    that the first category is zero for all items, treated null categories as
    specified via argument \code{"nullcats"} and without observations with zero
    weight,}
  \item{items}{logical vector of length \code{ncol(y)}, indicating
    which items were used during estimation,}
  \item{categories}{list of length \code{ncol(y)}, containing integer
    vectors starting from one to the number of categories minus one per
    item,}
  \item{n}{number of observations (with non-zero weights),}
  \item{n_org}{original number of observations in \code{y},}
  \item{weights}{the weights used (if any),}
  \item{na}{logical indicating whether the data contain \code{NA}s,}
  \item{nullcats}{currently always \code{NULL} as eventual items with null
    categories are handled via \code{"downcode"},}
  \item{impact}{either \code{NULL} or the supplied \code{impact} variable
    with the levels reordered in decreasing order (if this has not been the case
    prior to fitting the model),}
  \item{loglik}{log-likelihood of the fitted model,}
  \item{df}{number of estimated (more precisely, returned) model parameters,}
  \item{code}{convergence code from \code{mirt},}
  \item{iterations}{number of iterations used by \code{mirt},}
  \item{reltol}{convergence threshold passed to \code{mirt},}
  \item{grouppars}{the logical \code{grouppars} value,}
  \item{type}{the \code{type} of model restriction specified,}
  \item{mirt}{the \code{mirt} object fitted internally,}
  \item{call}{original function call.}
}

\references{
  Baker FB, Kim SH (2004).
    \emph{Item Response Theory: Parameter Estimation Techniques}.
    Chapman & Hall/CRC, Boca Raton.

  Bock RD, Aitkin M (1981).
    Marginal Maximum Likelihood Estimation of Item Parameters: Application of an EM Algorithm.
    \emph{Psychometrika}, \bold{46}(4), 443--459.

  Chalmers RP (2012).
    mirt: A Multidimensional Item Response Theory Package for the R Environment.
    \emph{Journal of Statistical Software}, \bold{48}(6), 1--29.
    \doi{10.18637/jss.v048.i06}

  Debelak R, Strobl C (2019).
    Investigating Measurement Invariance by Means of Parameter Instability Tests for 2PL and 3PL Models.
    \emph{Educational and Psychological Measurement}, \bold{79}(2), 385--398.
    \doi{10.1177/0013164418777784}

  Muraki E (1992).
    A Generalized Partial Credit Model: Application of an EM Algorithm.
    \emph{Applied Psychological Measurement}, \bold{16}(2), 159--176.

  Schneider L, Strobl C, Zeileis A, Debelak R (2022).
    An R Toolbox for Score-Based Measurement Invariance Tests in IRT Models.
    \emph{Behavior Research Methods}, forthcoming.
    \doi{10.3758/s13428-021-01689-0}
}

\seealso{\code{\link{pcmodel}}, \code{\link{rsmodel}}, \code{\link{nplmodel}},
  \code{\link{raschmodel}}, \code{\link{btmodel}}}

\examples{
if(requireNamespace("mirt")) {

o <- options(digits = 4)

## mathematics 101 exam results
data("MathExam14W", package = "psychotools")

## generalized partial credit model
gpcm <- gpcmodel(y = MathExam14W$credit)
summary(gpcm)

## how to specify starting values as a vector of model parameters
st <- coef(gpcm)
gpcm <- gpcmodel(y = MathExam14W$credit, start = st)
## or a list containing a vector of slopes and a list of intercept vectors
## itemwise
set.seed(0)
st <- list(a = rlnorm(13, 0, 0.0625), d = replicate(13, rnorm(2, 0, 1), FALSE))
gpcm <- gpcmodel(y = MathExam14W$credit, start = st)

## visualizations
plot(gpcm, type = "profile")
plot(gpcm, type = "regions")
plot(gpcm, type = "piplot")
plot(gpcm, type = "curves", xlim = c(-6, 6))
plot(gpcm, type = "information", xlim = c(-6, 6))
## visualizing the IRT parametrization
plot(gpcm, type = "curves", xlim = c(-6, 6), items = 1)
abline(v = threshpar(gpcm)[[1]])
abline(v = itempar(gpcm)[1], lty = 2)

options(digits = o$digits)
}
}

\keyword{regression}

