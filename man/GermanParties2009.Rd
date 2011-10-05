\name{GermanParties2009}
\alias{GermanParties2009}
\encoding{latin1}

\title{Choice among German Political Parties}

\description{
  Preferences of 192 respondents choosing among five German political parties
  and abstention from voting.
}

\usage{data("GermanParties2009")}

\format{
A data frame containing 192 observations on 6 variables.
\describe{
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    All 15 pairwise choices among five German parties and abstention from
    voting.}
  \item{ordered.pref}{Ordered paired comparison of class
    \code{\link{paircomp}}. Same as \code{preference}, but within-pair order
    is recognized.}
  \item{gender}{Factor coding gender.}
  \item{age}{Integer. Age of the respondents in years.}
  \item{education}{Ordered factor. Level of education: 1 no degree, 2
    Hauptschule (Secondary General School), 3 Realschule (Intermediate
    Secondary School), 4 Gymnasium (High School), 5 Studium (University)}
  \item{crisis}{Factor. Do you feel affected by the economic crisis?}
  \item{interviewer}{Factor. Interviewer id.}
}
}

\details{
  A survey was conducted at the Department of Psychology, Universität
  Tübingen, in June 2009, three months before the German election. The sample
  was stratified by gender and age (younger versus older than 30 years) with
  48 participants in each group.

  The parties to be compared were Die Linke (socialists), Die Grünen
  (ecologists), SPD (social democrats), CDU/CSU (conservatives), and FDP
  (liberals). In addition, there was the option of abstaining from voting
  (coded as \code{none}).

  Participants were presented with 15 pairs of options in random order. On
  each trial, their task was to choose the party they would rather vote for at
  an election for the German parliament. A given participant compared each
  pair in one order only, hence the NA's in \code{ordered.pref}.

  In order to minimize response biases, the pairs of options were read to the
  participants one at a time. Participants made their choices by crossing
  either \dQuote{First Option} or \dQuote{Second Option} on an anonymous
  response sheet.

  The interviewers were Psychology Master's students who collected the data
  for course credits. Since they mainly interviewed people they knew, the
  results are not representative of the political opinions in Germany. As far
  as the winner of the survey (Die Grünen) is concerned, however, the results
  agree with the outcome of the election for the Tübingen voters.

  The results of the election on September 27, 2009 (number of
  so-called Zweitstimmen in percent) were:

  \tabular{lrr}{
             \tab Germany \tab Tübingen\cr
  Die Linke  \tab    11.9 \tab      8.5\cr
  Die Grünen \tab    10.7 \tab     27.9\cr
  SPD        \tab    23.0 \tab     21.1\cr
  CDU/CSU    \tab    33.8 \tab     23.0\cr
  FDP        \tab    14.6 \tab     13.9\cr
  Others     \tab     6.0 \tab      5.7\cr
  }

  The voter turnout was 70.8 percent in Germany and 80.5 percent in Tübingen.
}

% \references{
%   FIX ME
% }

\seealso{\code{\link{paircomp}}}

\examples{
data("GermanParties2009", package = "psychotools")
summary(GermanParties2009$preference)
}

\keyword{datasets}
