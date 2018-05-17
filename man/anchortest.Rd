\name{anchortest}  
\alias{anchortest}
\alias{anchortest.default}
\alias{anchortest.formula}
\alias{print.anchortest}
\alias{print.summary.anchortest}
\alias{summary.anchortest}

\title{
  Anchor methods for the detection of uniform DIF in the Rasch model
}

\description{
  The \code{anchortest} function provides a Wald test (see,
  e.g., Glas, Verhelst, 1995) for the detection of uniform differential
  item functioning (DIF) in the Rasch model between two pre-specified
  groups. A variety of anchor methods is available to build a common
  scale necessary for the comparison of the item parameters in the Rasch
  model.
}

\usage{
anchortest(object, \dots)
\method{anchortest}{default}(object, object2,
  class = c("constant", "forward", "all-other", "fixed"), select = NULL,
  test = TRUE, adjust = "none", length = NULL, range = c(0.1, 0.8), \dots)
\method{anchortest}{formula}(formula, data = NULL, subset = NULL,
  na.action = NULL, weights = NULL, model = raschmodel, \dots)
}

\arguments{
  \item{object, object2}{Fitted model objects of class ``raschmodel''
    estimated via conditional maximum likelihood using \code{\link{raschmodel}}.}
  \item{\dots}{further arguments passed over to an internal call
    of \code{\link{anchor.default}} in the formula method. In the
  default method, these additional arguments are currently not being used.}
  \item{class}{character. Available anchor classes are the
    \code{constant} anchor class implying a constant anchor length
    defined by \code{length}, the iterative \code{forward} anchor class
    that iteratively includes items in the anchor and
    the \code{all-other} anchor class, for an overview see Kopf et al. (2015a).
    Additionally, the class can be \code{fixed}, then \code{select}
    needs to be the numeric index of the fixed selected anchor items.}
  \item{select}{character or numeric. Several anchor selection strategies are
    available: \code{"MTT"}, \code{"MPT"}, \code{"MT"}, \code{"MP"},
    \code{"NST"}, \code{"AO"}, \code{"AOP"}. For details see \code{\link{anchor}}.
    Alternatively, for \code{class = "fixed"}, \code{select} needs to be the
    numeric index of the fixed selected anchor items. Defaults are set such
    that \code{class = "constant"} is combined with \code{selection = "MPT"}
    while \code{class = "forward"} is combined with \code{selection = "MTT"}.
    And if \code{select} is numeric, then \code{class = "fixed"} is used.}
  \item{test}{logical. Should the Wald test be returned for the
    intended anchor method as final DIF test?}
  \item{adjust}{character. Should the final DIF test be adjusted for
    multiple testing? For the type of adjustment,
    see \code{\link[multcomp]{summary.glht}} and \code{\link{p.adjust}}.}
  \item{length}{integer. It pre-defines a maximum anchor length.
    Per default, the \code{forward} anchor grows up to the proportion
    of currently presumed DIF-free items specified in \code{range} and
    the \code{constant} anchor class selects four anchor items, unless
    an explicit limiting number is defined in \code{length} by the user.}
  \item{range}{numeric vector of length 2. The first element is
    the percentage of first anchor candidates to be excluded for
    consideration when the \code{forward} anchor class is used and the
    second element determines a percentage of currently presumed DIF-free
    items up to which the anchor from the \code{forward} anchor class is
    allowed to grow.}
    \item{formula}{formula of type \code{y ~ x} where \code{y} specifies a
    matrix of dichotomous item responses and \code{x} the grouping
    variable, e.g., gender, for which DIF should be tested for.}
  \item{data}{a data frame containing the variables of the specified
    \code{formula}.}
  \item{subset}{logical expression indicating elements or rows to keep:
    missing values are taken as false.}
  \item{na.action}{a function which indicates what should happen when the data
    contain missing values (\code{NA}s).}
  \item{weights}{an optional vector of weights (interpreted as case weights).}
  \item{model}{an IRT model fitting function with a suitable \code{itempar}
    method, by default \code{\link{raschmodel}}.}
}

\details{
  To conduct the Wald test (see, e.g., Glas, Verhelst, 1995) for uniform 
  DIF in the Rasch model, the user needs to specify an anchor method. 
  The anchor methods can be divided in an anchor class that determines 
  characteristics of the anchor method and an anchor selection that 
  determines the ranking order of candidate anchor items.

  Explicit anchor selection strategies are used in the \code{constant}
  anchor class and in the iterative \code{forward} anchor class, for a
  detailed description see \code{\link{anchor}}. Since \eqn{k-1}{k-1}
  parameters are free in the estimation, only \eqn{k-1}{k-1} estimated
  standard errors result. Thus, the first anchor item obtains no DIF test
  result and we report \eqn{k-1}{k-1} test results. This decision is
  applied only to those methods that rely on an explicit anchor
  selection strategy.

  In the \code{constant} anchor class, the anchor length is pre-defined
  by the user within the argument \code{length}. The default is four
  anchor items. The iterative \code{forward} class starts with a single
  anchor item and includes items in the anchor as long as the anchor
  length is shorter than a certain percentage of the number of items that
  do not display statistically significant DIF. The default proportion is
  set to 0.8 in the argument \code{range}. Alternatively, the user is
  allowed to set a maximum number of anchor items using the argument
  \code{length}. Both anchor classes require an explicit anchor selection
  strategy as opposed to the \code{all-other} anchor class.

  The \code{all-other} anchor class is here not considered as explicit
  anchor selection and, thus, only included in the \code{anchortest}
  function. For the \code{all-other} anchor class, the strategy is set to
  \code{"none"}, since all items except for the item currently studied
  for DIF are used as anchor. Thus, no explicit anchor selection strategy
  is required and we report \eqn{k}{k} test results. Note that the \code{all-other} 
  anchor class requires strong prior knowledge that DIF is balanced.

  See Kopf et al. (2015ab) for a detailed introduction. For convenience
  a trivial \code{"fixed"} anchor class is provided where the \code{select}ed
  anchor is given directly (e.g., as chosen by a practitioner or by some
  other anchor selection method).
}

\value{
  An object of class \code{anchor}, i.e. a list including
  \item{anchor_items}{the anchor items for DIF analysis.}
  \item{ranking_order}{a ranking order of candidate anchor items.}
  \item{criteria}{the criterion values obtained by the respective anchor
    selection.}
  \item{anchored_item_parameters}{ the anchored item parameters using
    the anchor items.}
  \item{anchored_covariances}{the anchored covariance matrices using
    the anchor items.}
  \item{final_tests}{the final Wald test for uniform DIF detection if
    intended.}
}

\references{
  Glas CAW, Verhelst ND (1995).
    \dQuote{Testing the Rasch Model.}
    In Fischer GH, Molenaar IW (eds.),
    \emph{Rasch Models: Foundations, Recent Developments, and Applications}, chapter 5.
    Springer-Verlag, New York.

  Kopf J, Zeileis A, Strobl C (2015a).
    A Framework for Anchor Methods and an Iterative Forward Approach for DIF Detection.
    \emph{Applied Psychological Measurement}, \bold{39}(2), 83--103.
    \doi{10.1177/0146621614544195}

  Kopf J, Zeileis A, Strobl C (2015b).
    Anchor Selection Strategies for DIF Analysis: Review, Assessment, and New Approaches.
    \emph{Educational and Psychological Measurement}, \bold{75}(1), 22--56.
    \doi{10.1177/0013164414529792}

  Wang WC (2004).
    Effects of Anchor Item Methods on the Detection of Differential Item Functioning within the Family of Rasch Models.
    \emph{Journal of Experimental Education}, \bold{72}(3), 221--261.

  Woods C (2009).
    Empirical Selection of Anchors for Tests of Differential Item Functioning.
    \emph{Applied Psychological Measurement}, \bold{33}(1), 42--57.
}

\seealso{\code{\link{anchor}}}

\examples{
if(requireNamespace("multcomp")) {

o <- options(digits = 4)

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Rasch model for the self-to-blame situations; gender DIF test
raschmodels <- with(VerbalAggression, lapply(levels(gender), function(i) 
  raschmodel(resp2[gender == i, 1:12])))

## four anchor items from constant anchor class using MPT-selection
const1 <- anchortest(object = raschmodels[[1]], object2 = raschmodels[[2]],
  class = "constant", select = "MPT", length = 4)
const1
summary(const1)

## iterative forward anchor class using MTT-selection
set.seed(1)
forw1 <- anchortest(object = raschmodels[[1]], object2 = raschmodels[[2]], 
  class = "forward", select = "MTT", test = TRUE,
  adjust = "none", range = c(0.05,1))
forw1

## the same using the formula interface
set.seed(1)
forw2 <- anchortest(resp2[, 1:12] ~ gender, data = VerbalAggression,
  class = "forward", select = "MTT", test = TRUE,
  adjust = "none", range = c(0.05, 1))

## compare both results
all.equal(forw1, forw2, check.attributes = FALSE)

## DIF test with fixed given anchor (arbitrarily selected to be items 1 and 2)
anchortest(object = raschmodels[[1]], object2 = raschmodels[[2]], select = 1:2)

options(digits = o$digits)
}
}

\keyword{regression}
