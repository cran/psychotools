\name{VerbalAggression}
\alias{VerbalAggression}

\title{Situation-Response Questionnaire on Verbal Aggression}

\description{
  Responses of 316 subjects to 24 items describing possible
  reactions to 4 different frustrating situations. 
}

\usage{data("VerbalAggression")}

\format{
A data frame containing 316 observations on 4 variables.
\describe{
  \item{resp}{Item response matrix with values 0/1/2 coding no/perhaps/yes, respectively.}
  \item{resp2}{Dichotomized item response matrix with perhaps/yes merged to 1.}
  \item{gender}{Factor coding gender.}
  \item{anger}{Trait anger, assessed by the Dutch adaptation
    of the state-trait anger scale (STAS).}
}
}

\details{
The 24 items are constructed by factorial combination of four different
frustrating situations (see below), three possible verbally aggressive responses
(curse, scold, shout), and two behavioural models (want, do). The four situations are
\tabular{ll}{
  S1:\tab A bus fails to stop for me.\cr
  S2:\tab I miss a train because a clerk gave me faulty information.\cr
  S3:\tab The grocery store closes just as I am about to enter.\cr
  S4:\tab The operator disconnects me when I used up my last 10 cents for a call.\cr
}
Note that the first two situations are other-to-blame situations, and the latter two are self-to-blame
situations.

The subjects were 316 first-year psychology students from a university in the Dutch
speaking part of Belgium. Participation was a partial fulfillment of the requirement
to participate in research. The sample consists of 73 males and 243 females,
reflecting the gender proportion among psychology students. The average
age was 18.4.
}

\source{
Online materials accompanying De Boeck and Wilson (2004).

\url{http://bear.soe.berkeley.edu/EIRM/}
}

\references{
De Boeck, P., Wilson, M. (eds) (2004).
  Explanatory Item Response Models: A Generalized Linear and Nonlinear Approach.
  New York: Springer-Verlag.

Smits, D.J.M., De Boeck, P., Vansteelandt, K. (2004).
  The Inhibition of Verbally Aggressive Behaviour
  \emph{European Journal of Personality}, \bold{18}, 537-555.
  doi:10.1002/per.529
}

\seealso{\code{\link{raschmodel}}}

\examples{
data("VerbalAggression", package = "psychotools")

## Rasch model for the self-to-blame situations
m <- raschmodel(VerbalAggression$resp2[, 1:12])
summary(m)
plot(m)
}

\keyword{datasets}
