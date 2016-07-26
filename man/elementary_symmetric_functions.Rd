\name{elementary_symmetric_functions}
\alias{elementary_symmetric_functions}

\title{Calculation of the Elementary Symmetric Functions and Their
  Derivatives}

\description{
  Calculation of \code{elementary_symmetric_functions} (ESFs), their first and,
  in the case of dichotomous items, second derivatives with sum or
  difference algorithm for the Rasch, rating scale and partial credit
  model.
}

\usage{
elementary_symmetric_functions(par, order = 0L, log = TRUE,
  diff = FALSE, engine = NULL)
}

\arguments{
  \item{par}{numeric vector or a list. Either a vector of item difficulty
    parameters of dichotomous items (Rasch  model) or a list of
    item-category parameters of polytomous items (rating scale and
    partial credit model).}
  \item{order}{integer between 0 and 2, specifying up to which derivative
    the ESFs should be calculated. Please note, second order derivatives
    are currently only possible for dichtomous items in an R
    implementation \code{engine == "R".}}
  \item{log}{logical. Are the parameters given in \code{par} on log
    scale? Primarily used for internal recursive calls of
    \code{elementary_symmetric_functions}.}  
  \item{diff}{logical. Should the first and second derivatives (if
    requested) of the ESFs calculated with sum (\code{FALSE})
    or difference algorithm (\code{TRUE}).}
  \item{engine}{character, either \code{"C"} or \code{"R"}. If the
    former, a C implementation is used to calculcate the ESFs and their
    derivatives, otherwise (\code{"R"}) pure R code is used.}  
}

\value{
  \code{elementary_symmetric_function} returns a list of length 1 + \code{order}.

  If \code{order = 0}, then the first (and only) element is a numeric
  vector with the ESFs of order 0 to the maximum score possible with
  the given parameters.

  If \code{order = 1}, the second element of the list contains a
  matrix, with the rows corresponding to the possible scores and the
  columns corresponding to the derivatives with respect to the i-th
  parameter of \code{par}.

  For dichotomous items and \code{order = 2}, the third element of the
  list contains an array with the second derivatives with respect to
  every possible combination of two parameters given in \code{par}. The
  rows of the individual matrices still correspond to the possibles
  scores (orders) starting from zero.
}

\details{
  Depending on the type of \code{par}, the elementary symmetric
  functions for dichotomous (\code{par} is a numeric vector) or
  polytomous items (\code{par} is a list) are calculated.
  
  For dichotomous items, the summation and difference algorithm
  published in Liou (1994) is used. For calculating the second order
  derivatives, the equations proposed by Jansens (1984) are employed.

  For polytomous items, the summation and difference algorithm published
  by Fischer and Pococny (1994) is used (see also Fischer and Pococny,
  1995).
}

\references{
  Liou M (1994).
    More on the Computation of Higher-Order Derivatives of the Elementary Symmetric Functions in the Rasch Model.
    \emph{Applied Psychological Measurement}, \bold{18}, 53--62.

  Jansen PGW (1984).
    Computing the Second-Order Derivatives of the Symmetric Functions in the Rasch Model.
    \emph{Kwantitatieve Methoden}, \bold{13}, 131--147.

  Fischer GH, and Ponocny I (1994).
    An Extension of the Partial Credit Model with an Application to the Measurement of Change.
    \emph{Psychometrika}, \bold{59}(2), 177--192.

  Fischer GH, and Ponocny I (1995).
    \dQuote{Extended Rating Scale and Partial Credit Models for Assessing Change.}
    In Fischer GH, and Molenaar IW (eds.).
    \emph{Rasch Models: Foundations, Recent Developments, and Applications.} 
}

\examples{
\donttest{
 ## zero and first order derivatives of 100 dichotomous items
 di <- rnorm(100)
 system.time(esfC <- elementary_symmetric_functions(di, order = 1))
 
 ## again with R implementation
 system.time(esfR <- elementary_symmetric_functions(di, order = 1,
 engine = "R"))

 ## are the results equal?
 all.equal(esfC, esfR)
}

 ## calculate zero and first order elementary symmetric functions
 ## for 10 polytomous items with three categories each.
 pi <- split(rnorm(20), rep(1:10, each = 2))
 x <- elementary_symmetric_functions(pi)

 ## use difference algorithm instead and compare results
 y <- elementary_symmetric_functions(pi, diff = TRUE)
 all.equal(x, y)
}

\keyword{misc}

