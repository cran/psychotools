\name{SourceMonitoring}
\alias{SourceMonitoring}

\title{Performance in a Source-Monitoring Experiment}

\description{
  Response frequencies of 128 participants who took part in a
  source-monitoring experiment with two sources.
}

\usage{data("SourceMonitoring")}

\format{
A data frame containing 128 observations on 15 variables.
\describe{
  \item{sources}{Factor. Sources A and B.}
  \item{age}{Integer. Age of the respondents in years.}
  \item{gender}{Factor coding gender.}
  \item{aa}{Number of source A items judged to be of source A.}
  \item{ab}{Number of source A items judged to be of source B.}
  \item{an}{Number of source A items judged to be new.}
  \item{ba}{Number of source B items judged to be of source A.}
  \item{bb}{Number of source B items judged to be of source B.}
  \item{bn}{Number of source B items judged to be new.}
  \item{na}{Number of new items judged to be of source A.}
  \item{nb}{Number of new items judged to be of source B.}
  \item{nn}{Number of new items judged to be new.}
}
}

\details{
  In a source-monitoring experiment with two sources, participants study
  items from two different sources, A and B. The final memory test consists
  of A and B items along with new distractor items, N. Participants are
  required to classify each item as A, B, or N.

  In an experiment at the Department of Psychology, University of Tuebingen,
  two source conditions were used in the study phase: In the first one,
  subjects had to read items either quietly (think) or aloud (say). In the
  second, they wrote them down (write) or read them aloud (say).

  The data were analyzed using the multinomial processing tree model of
  source monitoring (Batchelder & Riefer, 1990).
}

\source{
  Wickelmaier, F., and Zeileis, A. (2013). A first implementation of
  recursive partitioning for multinomial processing tree models. Presented
  at the Psychoco 2013 International Workshop on Psychometric Computing,
  February 14-15, Zurich, Switzerland. 
}

\references{
  Batchelder, W.H., and Riefer, D.M. (1990). Multinomial processing tree
  models of source monitoring. Psychological Review, 97, 548-564.
}

\examples{
data("SourceMonitoring", package = "psychotools")
xtabs(~ gender + I(age >= 30) + sources, SourceMonitoring)
}

\keyword{datasets}
