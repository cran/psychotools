\name{MemoryDeficits}
\alias{MemoryDeficits}

\title{Memory Deficits in Psychiatric Patients}

\description{
  Response frequencies of 96 patients who took part in a pair-clustering
  experiment to assess their memory deficits.
}

\usage{data("MemoryDeficits")}

\format{
A data frame containing 576 observations on 7 variables.
\describe{
  \item{ID}{Participant ID.}
  \item{group}{Factor with four levels specifying patient or control group
    of participant.}
  \item{trial}{Trial number from 1 to 6.}
  \item{E1}{Number of pairs recalled adjacently.}
  \item{E2}{Number of pairs recalled non-adjacently.}
  \item{E3}{Number of single pair members recalled.}
  \item{E4}{Number of non-recalled pairs.}
}
}

\details{
  Riefer, Knapp, Batchelder, Bamber and Manifold (2002) report a study on
  memory deficits in schizophrenic (n = 29) and organic alcoholic (n = 21)
  patients who were compared to two matched control groups (n = 25, n = 21).
  Participants were presented with 20 pairs of semantically related words. In
  a later memory test, they freely recalled the presented words. This
  procedure was repeated for a total of six study and test trials. Responses
  were classified into four categories: each pair is recalled adjacently (E1),
  each pair is recalled non-adjacently (E2), one word in a pair is recalled
  (E3), neither word in a pair is recalled (E4).
}

\source{
  The data were made available by William H. Batchelder.
}

\references{
  Riefer, D.M., Knapp, B.R., Batchelder, W.H., Bamber, D., & Manifold, V.
  (2002). Cognitive psychometrics: Assessing storage and retrieval deficits in
  special populations with multinomial processing tree models. Psychological
  Assessment, 14, 184-201.
}

\examples{
data("MemoryDeficits", package = "psychotools")
aggregate(cbind(E1, E2, E3, E4) ~ trial + group, MemoryDeficits, sum)
}

\keyword{datasets}
