\name{anchor}
\alias{anchor}
\alias{anchor.default}
\alias{anchor.formula}
\alias{print.anchor}
\alias{print.summary.anchor}
\alias{summary.anchor}
\encoding{latin1}

\title{Anchor Methods for the Detection of Uniform DIF the Rasch Model}

\description{
  The \code{anchor} function provides a variety of anchor
  methods for the detection of uniform differential item functioning (DIF)
  in the Rasch model between two pre-specified groups. These methods can
  be divided in an anchor class that determines characteristics of the
  anchor method and an anchor selection that determines the ranking order
  of candidate anchor items. The aim of the \code{anchor} function is to
  provide anchor items for DIF testing, e.g. with
  \code{\link{anchortest}}.
}

\usage{
anchor(object, \dots)
\method{anchor}{default}(object, object2,
  class = c("constant", "forward"), select = NULL,
  length = NULL, range = c(0.1, 0.8), \dots)
\method{anchor}{formula}(formula, data = NULL, subset = NULL,
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
    defined by \code{length} and the iterative \code{forward} anchor
    class, for an overview see Kopf et al. (2015a).}
  \item{select}{character. Several anchor selection strategies are
    available: \code{"MTT"}, \code{"MPT"}, \code{"MT"}, \code{"MP"},
    \code{"NST"}, \code{"AO"}, \code{"AOP"} based on
    Kopf et al. (2015b) as well as \code{"Gini"}, \code{"CLF"},
    \code{"GiniT"}, \code{"CLFT"} based on Strobl et al. (2021).
    The latter four can only be combined with \code{class = "constant"}
    and \code{length = 1}. The default is \code{select = "Gini"} unless
    either \code{length > 1} where \code{"MPT"} is used
    or \code{class = "forward"} where \code{"MTT"} is used. For more
    details see below.}
  \item{length}{integer. It pre-defines a maximum anchor length.
    Per default, the \code{forward} anchor grows up to the proportion
    of currently presumed DIF-free items specified in \code{range}
    and the \code{constant} anchor class selects one anchor item,
    unless an explicit limiting number is defined in \code{length}
    by the user.}
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
The anchor methods provided consist of an anchor class that determines
characteristics of the anchor and an anchor selection that determines
the ranking order of candidate anchor items.

In the \code{constant} anchor class, the anchor length is pre-defined by the
user within the argument \code{length}, defaulting to a length of one.
In contrast, the iterative \code{forward} class starts with a single anchor item and
includes items in the anchor as long as the anchor length is shorter
than a certain percentage of the number of items that do not display
statistically significant DIF (default: 0.8). Furthermore, a percentage 
of first anchor candidates is excluded from consideration (default: 0.1) 
and the user is allowed to set a maximum number of anchor items using the 
argument \code{length}. A detailed description of the anchor classes can be found
in Kopf et al. (2015a).

In more recent work Strobl et al. (2021) suggest a simpler yet powerful 
anchor method based on inequality criteria like the Gini coefficient. A
similar approach based on the component loss function (CLF) was suggested
by Muthén & Asparouhov (2014). These criteria can be shown to attain their
optimium for a single-anchor, thus correponding to a \code{constant} class of
\code{length} 1. Due to the simple structure in combination with good
empirical performance the Gini-based selection was made the default in
version 0.7-0 of the package.

Both anchor classes require an explicit anchor selection strategy (as opposed to
the \code{all-other} anchor class which is therefore not included in the
function \code{anchor}). The anchor selection strategy determines the ranking order of
candidate anchor items. In case of two groups, each item \eqn{j, j = 1,
\ldots, k}{j, j = 1, \dots, k} (where \eqn{k}{k} denotes the number of
items in the test) obtains a criterion value \eqn{c_j}{c_j} that is
defined by the anchor selection strategy. The ranking order is
determined by the rank of the criterion value
rank\eqn{(c_j)}{(c_j)}. 

The criterion values \eqn{c_j}{c_j} for item \eqn{j}{j} from the different
anchor selection strategies are provided in the following equations:
\eqn{d_j}{d_j} denotes the difference of the item parameters,
\eqn{t_j}{t_j} the corresponding test statistic, and
\eqn{p_j}{p_j} the resulting p-values. In all cases, the anchor items
are given in parentheses. Furthermore,
\eqn{\mathrm{Gini}(\cdot)}{Gini(.)} denotes the Gini inequality index,
\eqn{\mathrm{CLF}(\cdot)}{CLF(.)} the component loss function (sum of square root values),
\eqn{1(\cdot)}{1(.)} the indicator function,
\eqn{\lceil 0.5\cdot k \rceil}{<0.5k>} the empirical 50\% quantile,
and \eqn{A_\mathrm{purified}}{A_{purified}} the anchor after purification steps.
More detailed descriptions are available in Strobl et al. (2021) and
Kopf et al. (2015b).

Gini selection (of item parameter differences) by Strobl et al. (2021):
\deqn{ c_j^\mathrm{Gini} = - \mathrm{Gini} (\{ |d_1(j)|, \ldots, |d_k(j)| \}) }{%
 c_j^{Gini} = -Gini({|d_1(j)|, \ldots, |d_k(j)|}) }

GiniT selection (of test statistics) similar to Strobl et al. (2021):
\deqn{ c_j^\mathrm{GiniT} = - \mathrm{Gini} (\{ |t_1(j)|, \ldots, |t_k(j)| \}) }{%
 c_j^{GiniT} = -Gini({|t_1(j)|, \ldots, |t_k(j)|}) }

CLF selection (of item parameter differences) by Muthén & Asparouhov (2014):
\deqn{ c_j^\mathrm{CLF} = \mathrm{CLF} (\{ |d_1(j)|, \ldots, |d_k(j)| \}) }{%
 c_j^{CLF} = CLF({|d_1(j)|, \ldots, |d_k(j)|}) }

CLFT selection (of test statistics) similar to Muthén & Asparouhov (2014):
\deqn{ c_j^\mathrm{CLFT} = \mathrm{CLF} (\{ |t_1(j)|, \ldots, |t_k(j)| \}) }{%
 c_j^{CLFT} = CLF({|t_1(j)|, \ldots, |t_k(j)|}) }

All-other selection by Woods (2009), here abbreviated AO:
\deqn{ c_j^\mathrm{AO} = | t_j (\{1,\ldots,k\}\setminus j) |}{%
 c_j^{AO} = | t_j ({1,\dots,k}\\j) |}

All-other purified selection by Wang et al. (2012), here abbreviated AOP:
\deqn{ c_j^\mathrm{AOP} = | t_j ( A_\mathrm{purified} ) |}{%
 c_j^{AOP} = | t_j ( A_{purified} ) |}

Number of significant threshold selection based on Wang et al. (2004), here abbreviated NST:
\deqn{ c_j^\mathrm{NST} = \sum_{l \in \{1,\ldots,k\} \setminus j} 1 \left\{ p_j ( \{l\} )  \leq \alpha \right\}) |}{%
 c_j^{NST} = \sum_{l in {1,\dots,k} \\ j} 1 { p_j ({l}) \le \alpha }) |}

Mean test statistic selection by Shih et al. (2009), here abbreviated MT:
\deqn{ c_j^\mathrm{MT} = \frac{1}{k-1} \sum_{l \in \{1,\ldots,k\} \setminus j} \left| t_j ( \{l\}) \right| }{%
 c_j^{MT} = 1/(k-1) \sum_{l in {1,\dots,k} \\ j} | t_j ({l}) |}

Mean p-value selection by Kopf et al. (2015b), here abbreviated MP:
\deqn{ c_j^\mathrm{MP} = - \frac{1}{k-1} \sum_{l \in \{1,\ldots,k\} \setminus j} p_j ( \{l\} ) }{%
c_j^{MP} = - 1/(k-1) \sum_{l in {1,\dots,k} \\ j} p_j ({l})}

Mean test statistic threshold selection by Kopf et al. (2015b), here abbreviated MTT:
\deqn{ c_j^\mathrm{MTT} = 
 \sum_{l \in \{1,\ldots,k\} \setminus j} 
1 \left\{ \left|  t_j ( \{l\} ) \right| >  
 \left( \left| \frac{1}{k-1} \sum_{l \in \{ 1, \ldots, k \} \setminus j} 
t_j ( \{l\} ) \right|
\right)_{\left( \lceil 0.5\cdot k\rceil \right)}
  \right\} }{%
c_j^\mathrm{MTT} = 
 \sum_{l in {1,\dots,k} \\ j} 
1 { |  t_j ({l}) | > ( | 1/(k-1) \sum_{l \in { 1,\dots, k } \\ j} 
t_j ({l}) | )_{<0.5k>}}  
}

Mean p-value threshold selection by Kopf et al. (2015b), here abbreviated MPT:
\deqn{ c_j^\mathrm{MPT} = 
 -  \sum_{l \in \{1,\ldots,k\} \setminus j} 
1 \left\{ p_j ( \{l\} ) >  \left( 
\frac{1}{k-1} \sum_{l \in \{ 1, \ldots, k \} \setminus j} 
 p_j ( \{l\} )
 \right)_{ \left( \lceil 0.5\cdot k\rceil \right)}
  \right\} }{%
c_j^\mathrm{MPT} = 
 -  \sum_{l in {1,\dots,k} \\ j} 
1 { p_j ({l}) > ( 1/(k-1) \sum_{l in { 1, \dots, k } \\ j} 
 p_j ({l}))_{<0.5k>}
  }
}

Kopf et al. (2015b) recommend to combine the \code{class = "constant"} with
\code{select = "MPT"} and the \code{class = "forward"} with \code{select = "MTT"},
respectively.

The \code{all-other} anchor class (that assumes that DIF is balanced i.e. no 
group has an advantage in the test) is here
not considered as explicit anchor selection and, thus, not included
in the \code{anchor} function (but in the \code{\link{anchortest}}
function). Note that the \code{all-other} anchor class requires strong prior
knowledge that DIF is balanced.
}

\value{
  An object of class \code{anchor}, i.e. a list including 
  \item{anchor_items}{the integer index for the selected anchor items}
  \item{ranking_order}{a ranking order (integer index) of the candidate anchor items by their criterion values}
  \item{criteria}{the criterion values obtained in the anchor selection for each item (unsorted)}
}

\references{
  Kopf J, Zeileis A, Strobl C (2015a).
    A Framework for Anchor Methods and an Iterative Forward Approach for DIF Detection.
    \emph{Applied Psychological Measurement}, \bold{39}(2), 83--103.
    \doi{10.1177/0146621614544195}

  Kopf J, Zeileis A, Strobl C (2015b).
    Anchor Selection Strategies for DIF Analysis: Review, Assessment, and New Approaches.
    \emph{Educational and Psychological Measurement}, \bold{75}(1), 22--56.
    \doi{10.1177/0013164414529792}

  Muthén B, Asparouhov T (2014).
    IRT Studies of Many Groups: The Alignment Method.
    \emph{Frontiers in Psychology}, \bold{5}, 978.
    \doi{10.3389/fpsyg.2014.00978}

  Shih CL, Wang WC (2009).
    Differential Item Functioning Detection Using the Multiple Indicators, Multiple Causes Method with a Pure Short Anchor.
    \emph{Applied Psychological Measurement}, \bold{33}(3), 184--199.

  Strobl C, Kopf J, Kohler L, von Oertzen T, Zeileis A (2021).
    Anchor Point Selection: Scale Alignment Based on an Inequality Criterion.
    \emph{Applied Psychological Measurement}, \bold{45}(3), 214--230.
    \doi{10.1177/0146621621990743}

  Wang WC (2004).
    Effects of Anchor Item Methods on the Detection of Differential Item Functioning within the Family of Rasch Models.
    \emph{Journal of Experimental Education}, \bold{72}(3), 221--261.

  Wang WC, Shih CL, Sun GW (2012).
    The DIF-Free-then-DIF Strategy for the Assessment of Differential Item Functioning.
    \emph{Educational and Psychological Measurement}, \bold{72}(4), 687--708.

  Woods C (2009).
    Empirical Selection of Anchors for Tests of Differential Item Functioning.
    \emph{Applied Psychological Measurement}, \bold{33}(1), 42--57.
}

\seealso{
\code{\link{anchortest}}
}

\examples{
## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Gini anchor (Strobl et al. 2021) for gender DIF in the self-to-blame situations
anchor(resp2[, 1:12] ~ gender , data = VerbalAggression)

## alternatively: based on fitted raschmodel objects
raschmodels <- with(VerbalAggression, lapply(levels(gender), function(i) 
  raschmodel(resp2[gender == i, 1:12])))
anchor(raschmodels[[1]], raschmodels[[2]])

if(requireNamespace("multcomp")) {

## four anchor items from constant anchor class using MPT-selection (Kopf et al. 2015b)
anchor(object = raschmodels[[1]], object2 = raschmodels[[2]], 
  class = "constant", select = "MPT", length = 4)

## iterative forward anchor class using MTT-selection (Kopf et al. 2015b)
set.seed(1)
fanchor <- anchor(object = raschmodels[[1]], object2 = raschmodels[[2]],
  class = "forward", select = "MTT", range = c(0.05, 1))
fanchor

## the same using the formula interface
set.seed(1)
fanchor2 <- anchor(resp2[, 1:12] ~ gender , data = VerbalAggression,
  class = "forward", select = "MTT", range = c(0.05, 1))

## criteria really the same?
all.equal(fanchor$criteria, fanchor2$criteria, check.attributes = FALSE)
}
}

\keyword{regression}
