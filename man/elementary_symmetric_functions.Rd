\name{elementary_symmetric_functions}
\alias{elementary_symmetric_functions}

\title{Calculation of the Elementary Symmetric Functions and Their
  Derivatives}

\description{
  Calculation of \code{elementary_symmetric_functions} (ESFs), their first and
  second derivatives with sum or difference algorithm for the Rasch model.
}

\usage{
elementary_symmetric_functions(par, order = 0, log = TRUE, diff = FALSE)
}

\arguments{
  \item{par}{numeric vector of either item difficulty parameters.}
  \item{order}{integer between 0 and 2, specifying up to which derivative
    the ESFs should be calculated.}  
  \item{log}{logical. Are the parameters given in \code{par} on log
    scale? Mainly used for internal recursive calls of
    \code{elementary_symmetric_functions}.}  
  \item{diff}{logical. Should the first derivatives (if requested) of
    the ESFs calculated with sum algorithm (\code{FALSE}) or difference
    algorithm (\code{TRUE}).}
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

  For \code{order = 2}, the third element of the
  list contains an array with the second derivatives with respect to
  every possible combination of two parameters given in \code{par}. The
  rows of the individual matrices still correspond to the possibles
  scores (orders) starting from zero.
}

\keyword{misc}
