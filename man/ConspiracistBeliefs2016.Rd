\name{ConspiracistBeliefs2016}
\alias{ConspiracistBeliefs2016}

\title{Generic Conspiracist Beliefs Scale (2016 Data)}

\description{
  Responses of 2449 persons to 15 five-point likert-rated items (0 = disagree to
  4 = agree) measuring belief in conspiracy theories as well as responses on 2
  covariates.
}

\usage{data("ConspiracistBeliefs2016", package = "psychotools")}

\format{
A data frame containing 2449 observations on 3 variables.
\describe{
  \item{resp}{Item response matrix with 15 items (see details below).}
  \item{area}{Factor coding the area one lived in as a child (\code{"rural"},
    \code{"suburban"}, \code{"urban"}).}
  \item{gender}{Factor coding gender (\code{"male"}, \code{"female"},
    \code{"other"}).}
}
}

\details{
The Open Source Psychometrics Project published this dataset collected online
in 2016. Persons responded to the Generic Conspiracist Beliefs (GCB) Scale
(Brotherton, French & Pickering, 2013) as well as other additional questions
primarily for personal amusement. At the end of the test but before the
results were displayed, users were asked if they would allow their responses to
be saved for research. Only users who agreed are part of this dataset.
Individuals with age lower than 13 years were not recorded. Moreover, two
persons stating their age to be 5555 years or higher as well as 44 persons
with missing data in \code{area} or \code{gender} were excluded from this
dataset. The 15 items of the GCB Scale are:
\tabular{rl}{
  Q1:\tab The government is involved in the murder of innocent citizens and/or
    well-known public figures, and keeps this a secret.\cr
  Q2:\tab The power held by heads of state is second to that of small unknown
    groups who really control world politics.\cr
  Q3:\tab Secret organizations communicate with extraterrestrials, but keep
    this fact from the public.\cr
  Q4:\tab The spread of certain viruses and/or diseases is the result of the
    deliberate, concealed efforts of some organization.\cr
  Q5:\tab Groups of scientists manipulate, fabricate, or suppress evidence in
    order to deceive the public.\cr
  Q6:\tab The government permits or perpetrates acts of terrorism on its own
    soil, disguising its involvement.\cr
  Q7:\tab A small, secret group of people is responsible for making all major
    world decisions, such as going to war.\cr
  Q8:\tab Evidence of alien contact is being concealed from the public.\cr
  Q9:\tab Technology with mind-control capacities is used on people without
    their knowledge.\cr
  Q10:\tab New and advanced technology which would harm current industry is
    being suppressed.\cr
  Q11:\tab  The government uses people as patsies to hide its involvement in
    criminal activity.\cr
  Q12:\tab Certain significant events have been the result of the activity of a
    small group who secretly manipulate world events.\cr
  Q13:\tab Some UFO sightings and rumors are planned or staged in order to
    distract the public from real alien contact.\cr
  Q14:\tab Experiments involving new drugs or technologies are routinely
    carried out on the public without their knowledge or consent.\cr
  Q15:\tab A lot of important information is deliberately concealed from the
    public out of self-interest.
}
Additional information can be found online (see below) via inspecting the
codebook contained in \file{GCBS.zip}.
}

\source{\url{https://openpsychometrics.org/_rawdata/}.}

\references{
  Brotherton R, French CC, Pickering AD (2013).
    Measuring Belief in Conspiracy Theories: The Generic Conspiracist Beliefs
    Scale.
    \emph{Frontiers in Psychology}, \bold{4}, 279.

  Open Source Psychometrics Project (2016).
    Data From: The Generic Conspiracist Beliefs Scale [Dataset].
    Retrieved from \url{https://openpsychometrics.org/_rawdata/}.
}

\seealso{\code{\link{gpcmodel}}}

\examples{
## overview
data("ConspiracistBeliefs2016", package = "psychotools")
str(ConspiracistBeliefs2016)

## response
plot(itemresp(ConspiracistBeliefs2016$resp))
## covariates
summary(ConspiracistBeliefs2016[, -1])
}

\keyword{datasets}

