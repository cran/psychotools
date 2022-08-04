\name{nplmodel}
\alias{nplmodel}
\alias{plmodel}

\alias{print.nplmodel}
\alias{summary.nplmodel}
\alias{print.summary.nplmodel}

\alias{coef.nplmodel}
\alias{confint.nplmodel}
\alias{bread.nplmodel}
\alias{estfun.nplmodel}
\alias{logLik.nplmodel}
\alias{vcov.nplmodel}

\title{Parametric Logistic Model (n-PL) Fitting Function}

\description{
  \code{nplmodel} is a basic fitting function for n-PL type parametric logistic IRT models
  (2PL, 3PL, 3PLu, 4PL, Rasch/1PL), providing a wrapper around
  \code{\link[mirt]{mirt}} and \code{\link[mirt]{multipleGroup}} relying on
  marginal maximum likelihood (MML) estimation via the standard EM algorithm.
}

\usage{
nplmodel(y, weights = NULL, impact = NULL,
  type = c("2PL", "3PL", "3PLu", "4PL", "1PL", "RM"),
  grouppars = FALSE, vcov = TRUE, 
  start = NULL, method = "BFGS", maxit = 500, reltol = 1e-5, \dots)
}

\arguments{
  \item{y}{item response object that can be coerced (via \code{\link[base]{as.matrix}})
    to a numeric matrix with scores 0, 1. Typically, either
    already a matrix, data frame, or dedicated object of class
    \code{\link{itemresp}}.}
  \item{weights}{an optional vector of weights (interpreted as case
    weights).}
  \item{impact}{an optional \code{factor} allowing for grouping
    the subjects (rows). If specified, a multiple-group model is fitted
    to account for impact (see details below). By default, no impact is
    modelled, i.e., a single-group model is used.}
  \item{type}{character string, specifying the type of parametric logistic
    IRT model to be estimated (see details below).}
  \item{grouppars}{logical. Should the estimated distributional group parameters
    of a multiple group model be included in the model parameters?}
  \item{vcov}{logical or character specifying the type of variance-covariance
    matrix (if any) computed for the final model. The default \code{vcov = TRUE}
    corresponds to \code{vcov = "Oakes"}, see \code{\link[mirt]{mirt}} for
    further options. If set to \code{vcov = FALSE} (or \code{vcov = "none"}),
    \code{vcov()} will return a matrix of \code{NA}s only.}
  \item{start}{an optional vector or list of starting values (see examples
    below).}
  \item{method, maxit, reltol}{control parameters for the optimizer employed
    by \code{\link[mirt]{mirt}} for the EM algorithm.}
  \item{\dots}{further arguments passed to \code{\link[mirt]{mirt}} or
    \code{\link[mirt]{multipleGroup}}, respectively.}
}

\details{
  \code{nplmodel} (\code{plmodel} for backward compatibility with earlier
  \pkg{psychotools} versions) provides a basic fitting function for n-PL type parametric logistic IRT
  models (2PL, 3PL, 3PLu, 4PL, Rasch/1PL) providing a wrapper around
  \code{\link[mirt]{mirt}} and \code{\link[mirt]{multipleGroup}} relying on
  MML estimation via the standard EM algorithm (Bock & Aitkin, 1981). Models are
  estimated under the slope/intercept parametrization, see e.g. Chalmers (2012).
  The probability of person \eqn{i} \sQuote{solving} item \eqn{j} is modelled as:
  \deqn{P(X_{ij} = 1|\theta_{i},a_{j},d_{j},g_{j},u_{j}) =
  g_{j} + \frac{(u_{j} - g_{j})}{1 + \exp{(-(a_{j}\theta_{i} + d_{j}))}}}

  A reparametrization of the intercepts to the classical IRT parametrization,
  \eqn{b_{j} = -\frac{d_{j}}{a_{j}}}, is provided via the corresponding
  \code{\link{itempar}} method.

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

  By default, \code{type} is set to \code{"2PL"}. Therefore, all so-called
  guessing parameters are fixed at 0 and all upper asymptotes are fixed at 1.
  \code{"3PL"} results in all upper asymptotes being fixed at 1 and \code{"3PLu"}
  results in all all guessing parameters being fixed at 0. \code{"4PL"} results
  in a full estimated model as specified above. Finally, if \code{type} is set to
  \code{"1PL"} (or equivalently \code{"RM"}), an MML-estimated Rasch model is
  being fitted. This means that all slopes are restricted to be equal across all
  items, all guessing parameters are fixed at 0 and all upper asymptotes are
  fixed at 1.

  Note that internally, the so-called guessing parameters and upper asymptotes
  are estimated on the logit scale (see also \code{\link[mirt]{mirt}}).
  Therefore, most of the basic methods below include a \code{logit} argument,
  which can be set to \code{TRUE} or \code{FALSE} allowing for a retransformation
  of the estimates and their variance-covariance matrix (if requested) using the
  logistic function and the delta method if \code{logit = FALSE}.

  \code{nplmodel} returns an object of class \code{"nplmodel"} for which
  several basic methods are available, including \code{print}, \code{plot},
  \code{summary}, \code{coef}, \code{vcov}, \code{logLik}, \code{estfun},
  \code{\link{discrpar}}, \code{\link{itempar}}, \code{\link{threshpar}},
  \code{\link{guesspar}}, \code{\link{upperpar}}, and \code{\link{personpar}}.
  
  Finally, if \code{type} is set to \code{"1PL"}, a Rasch model is estimated. 
  Here, a common slope parameter is estimated for all items, whereas the
  person parameters are assumed to follow a standard normal distribution.
  Please note that this variant of the Rasch model differs from the one used 
  by \code{\link[mirt]{mirt}}, which sets all slope parameters to 1, and 
  estimates the variance of the person parameters instead. Both variants 
  are mathematically equivalent and hence should lead to the same intercept parameter 
  estimates. For numerical reasons, \code{nplmodel} and \code{mirt} 
  can lead to slightly different item parameter estimates, though, under their 
  respective default settings, in particular when some items are very easy 
  or very difficult and the common slope parameter is large. A distinct advantage 
  of the variant used by \code{nplmodel} is that it allows a direct 
  comparison of the slope and intercept parameters with that estimated in more complex 
  IRT models, such as the 2PL model. 
}

\value{
  \code{nplmodel} returns an S3 object of class \code{"nplmodel"}, 
  i.e., a list of the following components:
  \item{coefficients}{estimated model parameters in slope/intercept parametrization,}
  \item{vcov}{covariance matrix of the model parameters,}
  \item{data}{modified data, used for model-fitting, i.e., without
    observations with zero weight,}
  \item{items}{logical vector of length \code{ncol(y)}, indicating
    which items were used during estimation,}
  \item{n}{number of observations (with non-zero weights),}
  \item{n_org}{original number of observations in \code{y},}
  \item{weights}{the weights used (if any),}
  \item{na}{logical indicating whether the data contain \code{NA}s,}
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
  \item{mirt}{the \code{mirt} object fitted internally.}
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

  Schneider L, Strobl C, Zeileis A, Debelak R (2022).
    An R Toolbox for Score-Based Measurement Invariance Tests in IRT Models.
    \emph{Behavior Research Methods}, forthcoming.
    \doi{10.3758/s13428-021-01689-0}
}

\seealso{\code{\link{raschmodel}}, \code{\link{gpcmodel}},
  \code{\link{rsmodel}}, \code{\link{pcmodel}}, \code{\link{btmodel}}}

\examples{
if(requireNamespace("mirt")) {

o <- options(digits = 4)

## mathematics 101 exam results
data("MathExam14W", package = "psychotools")

## 2PL
twopl <- nplmodel(y = MathExam14W$solved)
summary(twopl)

## how to specify starting values as a vector of model parameters
st <- coef(twopl)
twopl <- nplmodel(y = MathExam14W$solved, start = st)
## or a list containing a vector of slopes and a vector of intercepts
set.seed(0)
st <- list(a = rlnorm(13, 0, 0.0625), d = rnorm(13, 0, 1))
twopl <- nplmodel(y = MathExam14W$solved, start = st)

## visualizations
plot(twopl, type = "profile")
plot(twopl, type = "regions")
plot(twopl, type = "piplot")
plot(twopl, type = "curves", xlim = c(-6, 6))
plot(twopl, type = "information", xlim = c(-6, 6))
## visualizing the IRT parametrization
plot(twopl, type = "curves", xlim = c(-6, 6), items = 1)
abline(v = itempar(twopl)[1])
abline(h = 0.5, lty = 2)

## 2PL accounting for gender impact
table(MathExam14W$gender)
mtwopl <- nplmodel(y = MathExam14W$solved, impact = MathExam14W$gender,
  grouppars = TRUE)
summary(mtwopl)
plot(mtwopl, type = "piplot")
## specifying starting values as a vector of model parameters, note that in
## this example impact is being modelled and therefore grouppars must be TRUE
## to get all model parameters
st <- coef(mtwopl)
mtwopl <- nplmodel(y = MathExam14W$solved, impact = MathExam14W$gender,
  start = st)
## or a list containing a vector of slopes, a vector of intercepts and a vector
## of means and a vector of variances as the distributional group parameters
set.seed(1)
st <- list(a = rlnorm(13, 0, 0.0625), d = rnorm(13, 0, 1), m = 0, v = 1)
mtwopl <- nplmodel(y = MathExam14W$solved, impact = MathExam14W$gender,
  start = st)

## MML estimated Rasch model (1PL)
rm <- nplmodel(y = MathExam14W$solved, type = "1PL")
summary(rm)

options(digits = o$digits)
}
}

\keyword{regression}

