\name{PairClustering}
\alias{PairClustering}

\title{Pair Clustering Data in Klauer (2006)}

\description{
  Response frequencies of 63 participants who took part in a pair-clustering
  experiment.
}

\usage{data("PairClustering")}

\format{
A data frame containing 126 observations on 8 variables.
\describe{
  \item{ID}{Participant ID.}
  \item{trial}{Trial number, 1 or 2.}
  \item{E1}{Number of pairs recalled adjacently.}
  \item{E2}{Number of pairs recalled non-adjacently.}
  \item{E3}{Number of single pair members recalled.}
  \item{E4}{Number of non-recalled pairs.}
  \item{F1}{Number of recalled singleton words.}
  \item{F2}{Number of non-recalled singleton words.}
}
}

\details{
  Klauer (2006) reports a pair-clustering experiment with 63 participants,
  who were presented with ten pairs of related words and five unrelated
  singleton words. In a later memory test, they freely recalled the presented
  words. This procedure was repeated for two study and test trials. For pairs,
  responses were classified into four categories: both words in a pair are
  recalled adjacently (E1) or non-adjacently (E2), one word in a pair is
  recalled (E3), neither word in a pair is recalled (E4); for singletons, into
  two categories: word recalled (F1), word not recalled (F2).
}

\source{
  Stahl C, Klauer KC (2007).
    HMMTree: A Computer Program for Latent-Class Hierarchical Multinomial Processing Tree Models.
  \emph{Behavior Research Methods}, \bold{39}, 267--273.
}

\references{
  Klauer KC (2006).
    Hierarchical Multinomial Processing Tree Models: A Latent-Class Approach.
    \emph{Psychometrika}, \bold{71}, 1--31.
}

\examples{
data("PairClustering", package = "psychotools")
aggregate(cbind(E1, E2, E3, E4, F1, F2) ~ trial, PairClustering, sum)
}

\keyword{datasets}
