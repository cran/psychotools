\name{SourceMonitoring}
\alias{SourceMonitoring}

\title{Performance in a Source-Monitoring Experiment}

\description{
  Response frequencies of 128 participants who took part in a
  source-monitoring experiment with two sources.
}

\usage{data("SourceMonitoring")}

\format{
A data frame containing 128 observations on four components.
\describe{
  \item{sources}{Factor. Sources A and B.}
  \item{age}{Integer. Age of the respondents in years.}
  \item{gender}{Factor coding gender.}
  \item{y}{Matrix containing the response frequencies. The column names
    indicate the nine response categories:
    \tabular{rl}{
      \code{a.a} \tab Number of source A items judged to be of source A.\cr
      \code{a.b} \tab Number of source A items judged to be of source B.\cr
      \code{a.n} \tab Number of source A items judged to be new.\cr
      \code{b.a} \tab Number of source B items judged to be of source A.\cr
      \code{b.b} \tab Number of source B items judged to be of source B.\cr
      \code{b.n} \tab Number of source B items judged to be new.\cr
      \code{n.a} \tab Number of new items judged to be of source A.\cr
      \code{n.b} \tab Number of new items judged to be of source B.\cr
      \code{n.n} \tab Number of new items judged to be new.
    }
  }
}
}
\details{
  In a source-monitoring experiment with two sources, participants study
  items from two different sources, A and B. The final memory test consists
  of A and B items along with new distractor items, N. Participants are
  required to classify each item as A, B, or N.

  In an experiment at the Department of Psychology, University of Tuebingen
  (Wickelmaier & Zeileis, 2013, 2018), two source conditions were used in
  the study phase: Half of the subjects had to read items either quietly
  (source A = think) or aloud (source B = say). The other half had to write
  items down (source A = write) or read them aloud (source B = say).

  The data were analyzed using the multinomial processing tree model of
  source monitoring (Batchelder & Riefer, 1990).
}

\source{
  Wickelmaier F, Zeileis A (2013).
    A First Implementation of Recursive Partitioning for Multinomial
    Processing Tree Models.
    Presented at the Psychoco 2013 International Workshop on Psychometric
    Computing, February 14--15, Zurich, Switzerland. 
}

\references{
  Batchelder WH, Riefer DM (1990).
    Multinomial Processing Tree Models of Source Monitoring.
    \emph{Psychological Review}, \bold{97}, 548--564.

  Wickelmaier F, Zeileis A (2018).
    Using Recursive Partitioning to Account for Parameter Heterogeneity in
    Multinomial Processing Tree Models.
    \emph{Behavior Research Methods}, \bold{50}(3), 1217--1233.
    \doi{10.3758/s13428-017-0937-z}
}

\examples{
data("SourceMonitoring", package = "psychotools")
xtabs(~ gender + I(age >= 30) + sources, SourceMonitoring)
}

\keyword{datasets}
