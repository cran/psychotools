\name{Sim3PL}
\alias{Sim3PL}

\title{Simulated Data for fitting a 3PL and 3PLu}

\description{
  Simulated responses of 10000 persons to 10 dichotomous items under two
  different simulation conditions.
}

\usage{data("Sim3PL", package = "psychotools")}

\format{
A data frame containing 10000 observations on 2 variables.
\describe{
  \item{resp}{Item response matrix with 10 items (see details below).}
  \item{resp2}{Item response matrix with 10 items (see details below).}
}
}

\details{
  Data were simulated under the 3PL (\code{resp}) and 3PLu (\code{resp2}) (see
  \code{\link{nplmodel}}). For the 3PL scenario, the random number generator's
  seed was set to 277. For the 3PLu scenario, the random number generator's seed
  was set to 167. Person parameters \eqn{\theta_{i}} of 10000 persons were drawn
  from the standard normal distribution. Item difficulties \eqn{b_{j}} of 10
  items (under the classical IRT parametrization) were drawn from the standard
  normal distribution. Item discrimination parameters \eqn{a_{j}} were drawn
  from a log-normal distribution with a mean of \eqn{0} and a variance of
  \eqn{0.0625} on the log scale. For the 3PL, guessing parameters
  \eqn{g_{j}} were drawn from a uniform distribution with a lower limit of
  \eqn{0.1} and an upper limit of \eqn{0.2}. For the 3PLu, upper asymptote
  parameters \eqn{u_{j}} were drawn from a uniform distribution with a lower
  limit of \eqn{0.8} and an upper limit of \eqn{0.9}. In both scenarios, a
  \eqn{10000} x \eqn{10} matrix based on realizations of a uniform distribution
  with a lower limit of \eqn{0} and an upper limit of \eqn{1} was generated and
  compared to a \eqn{10000} x \eqn{10} matrix based on the probability function
  under the respective model. If the probability of person \eqn{i} solving item
  \eqn{j} exceeded the corresponding realization of the uniform distribution,
  this cell of the matrix was set to \eqn{1}, e.g., person \eqn{i} solved item
  \eqn{j}.
}

\seealso{\code{\link{nplmodel}}}

\examples{
## overview
data("Sim3PL", package = "psychotools")
str(Sim3PL)

## data generation
M <- 10000
N <- 10

## 3PL scenario
set.seed(277)
theta <- rnorm(M, 0, 1)
a <- rlnorm(N, 0, 0.25)
b <- rnorm(N, 0, 1)
g <- runif(N, 0.1, 0.2)
u <- rep(1, N)
probs <- matrix(g, M, N, byrow = TRUE) + matrix(u - g, M, N, byrow = TRUE) *
  plogis(matrix(a, M, N, byrow = TRUE) * outer(theta, b, "-"))
resp <- (probs > matrix(runif(M * N, 0, 1), M, N)) + 0
all.equal(resp, Sim3PL$resp, check.attributes = FALSE)

## 3PLu scenario
set.seed(167)
theta <- rnorm(M, 0, 1)
a <- rlnorm(N, 0, 0.25)
b <- rnorm(N, 0, 1)
g <- rep(0, N)
u <- runif(N, 0.8, 0.9)
probs <- matrix(g, M, N, byrow = TRUE) + matrix(u - g, M, N, byrow = TRUE) *
  plogis(matrix(a, M, N, byrow = TRUE) * outer(theta, b, "-"))
resp2 <- (probs > matrix(runif(M * N, 0, 1), M, N)) + 0
all.equal(resp2, Sim3PL$resp2, check.attributes = FALSE)
}

\keyword{datasets}

