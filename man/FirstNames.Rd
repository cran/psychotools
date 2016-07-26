\name{FirstNames}
\alias{FirstNames}
\encoding{latin1}

\title{Popularity of First Names}

\description{
  Preferences of 192 respondents choosing among six boys names with respect to
  their popularity.
}

\usage{data("FirstNames")}

\format{
A data frame containing 192 observations on 11 variables.
\describe{
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    All 15 pairwise choices among six boys names: Tim, Lucas, Michael, Robin,
    Benedikt, and Julius.}
  \item{ordered.pref}{Ordered paired comparison of class
    \code{\link{paircomp}}. Same as \code{preference}, but within-pair order
    is recognized.}
  \item{gender}{Factor coding gender.}
  \item{age}{Integer. Age of the respondents in years.}
  \item{education}{Ordered factor. Level of education: 1
    Hauptschule with degree (Secondary General School), 2 and 3 Realschule
    without and with degree (Intermediate Secondary School), 4 and 5
    Gymnasium without and with degree (High School), 6 and 7 Studium without
    and with degree (University).}
  \item{children}{Integer. Number of children.}
  \item{state}{Factor. State of Germany where participant grew up.}
  \item{state.reg}{Factor. The region (south, north-west, east) each state
    belongs to.}
  \item{fname}{Factor. Participant's fist name(s). (Umlaute in Jörg and Jürgen
    have been transliterated to Joerg and Juergen for portability of the data.)}
  \item{interviewer}{Factor. Interviewer id.}
  \item{gender.int}{Factor coding interviewer's gender.}
}
}

\details{
  A survey was conducted at the Department of Psychology, Universität
  Tübingen, in June 2009. The sample was stratified by gender and age (younger
  versus older than 30 years) with 48 participants in each group.  The
  interviewers were Psychology Master's students who collected the data for
  course credits.

  Participants were presented with 15 pairs of boys names in random order. On
  each trial, their task was to choose the name they would rather give to
  their own child.  The pairs of boys names were read to the participants one
  at a time.  A given participant compared each pair in one order only, hence
  the NA's in \code{ordered.pref}.

  The names were selected to fall within the upper (Tim, Lucas), mid (Michael,
  Robin) and lower (Benedikt, Julius) range of the top 100 of the most popular
  boys names in Germany in the years from 1990 to 1999
  (\url{http://www.beliebte-vornamen.de/3778-1990er-jahre.htm}). The names
  have either front (e, i) or back (o, u) vowels in the stressed syllables.
  Phonology of the name and attractiveness of a person have been shown to be
  related (Perfors, 2004; Hartung et al., 2009).
}

\references{
  Hartung F, Klenovsak D, Santiago dos Santos L, Strobl C, Zaefferer D (2009).
    Are Tims Hot and Toms Not? Probing the Effect of Sound Symbolism on Perception of Facial Attractiveness.
    Presented at the \emph{31th Annual Meeting of the Cognitive Science Society},
    July 27--August 1, Amsterdam, The Netherlands.
 
  Perfors A (2004).
    What's in a Name? The Effect of Sound Symbolism on Perception of Facial Attractiveness.
    Presented at the \emph{26th Annual Meeting of the Cognitive Science Society},
    August 5--7, Chicago, USA.
}

\seealso{\code{\link{paircomp}}}

\examples{
data("FirstNames", package = "psychotools")
summary(FirstNames$preference)
covariates(FirstNames$preference)
}

\keyword{datasets}
