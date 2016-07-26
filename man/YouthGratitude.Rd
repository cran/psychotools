\name{YouthGratitude}
\alias{YouthGratitude}

\title{Measuring Gratitude in Youth}

\description{
  Cross-section data on several gratitude scales for children and adolescents.
}

\usage{data("YouthGratitude")}

\format{
A data frame containing 1405 observations on 28 variables.
\describe{
  \item{id}{Integer person ID.}
  \item{age}{Age in years (10--19 years).}
  \item{agegroup}{Factor coding of age with levels \code{"10-11"}, \code{"12-13"},
    \code{"14"}, \code{"15"}, \code{"16"}, \code{"17-19"}.}

  \item{losd_1}{Life has been good to me.}
  \item{losd_2}{There never seems to be enough to go around, and I never seem to get my share. (Reverse scored.)}
  \item{losd_3}{I really don't think that I've gotten all the good things that I deserve in life. (Reverse scored.)}
  \item{losd_4}{More bad things have happened to me in my life than I deserve. (Reverse scored.)}
  \item{losd_5}{Because of what I've gone through in my life, I really feel like the world owes me something. (Reverse scored.)}
  \item{losd_6}{For some reason I never seem to get the advantages that others get. (Reverse scored.)}

  \item{sa_1}{Oftentimes I have been overwhelmed at the beauty of nature.}
  \item{sa_2}{Every Fall I really enjoy watching the leaves change colors.}
  \item{sa_3}{I think that it's important to 'Stop and smell the roses.'}
  \item{sa_4}{I think that it's important to pause often to 'count my blessings.'}
  \item{sa_5}{I think it's important to enjoy the simple things in life.}
  \item{sa_6}{I think it's important to appreciate each day that you are alive.}

  \item{ao_1}{I couldn't have gotten where I am today without the help of many people.}
  \item{ao_2}{Although I think it's important to feel good about your accomplishments, I think that it's also important to remember how others have contributed to my accomplishments.}
  \item{ao_3}{Although I'm basically in control of my life, I can't help but think about all those who have supported me and helped me along the way.}
  \item{ao_4}{I feel deeply appreciative for the things others have done for me in my life.}

  \item{gq6_1}{I have so much in life to be thankful for.}
  \item{gq6_2}{If I had to list everything that I felt thankful for,
    it would be a very long list.}
  \item{gq6_3}{When I look at the world, I don't see much to be thankful for.}
  \item{gq6_4}{I am thankful to a wide variety of people. (Reverse scored.)}
  \item{gq6_5}{As I get older I find myself more able to appreciate the people,
    events, and situations that have been part of my life history.}
  \item{gq6_6}{Long amounts of time can go by before I feel gratitude to something or someone. (Reverse scored.)}

  \item{gac_1}{Grateful.}
  \item{gac_2}{Thankful.}
  \item{gac_3}{Appreciative.}
}
}

\details{
The gratitude scales employed are:\cr
  GRAT: Gratitude, Resentment, Appreciation Test (1--9).\cr
    Short form with subscales LOSD (lack of a sense of deprivation),
    SA (simple appreciation), and AO (appreciation for others).\cr
  GQ-6: Gratitude Questionnaire-6 (1--7).\cr
  GAC: Gratitude Adjective Checklist (1--5).

The item \code{losd_1} has been omitted from all analyses in Froh et al. (2011)
because it loaded lowly on all factors. Hence \code{losd_1} is not listed in
Table B1 of Froh et al. (2011). Instead, the remaining items are labeled
\code{losd_1} to \code{losd_5}.
}

\source{
Provided by Jeff Froh and Jinyan Fan.
}

\references{
  Froh JJ, Fan J, Emmons RA, Bono G, Huebner ES, Watkins P (2011).
    Measuring Gratitude in Youth: Assessing the Psychometric Properties of Adult Gratitude Scales in Children and Adolescents.
    \emph{Psychological Assessment}, \bold{23}(2), 311--324.
}

\examples{
data("YouthGratitude", package = "psychotools")
summary(YouthGratitude)

## modeling can be carried out using package lavaan
\dontrun{
## remove cases with 'imputed' values (not in 1, ..., 9)
yg <- YouthGratitude[apply(YouthGratitude[, 4:28], 1, function(x) all(x %in% 1:9)), ]

## GQ-6
gq6_congeneric <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, group = "agegroup", meanstructure = TRUE)
gq6_tauequivalent <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = "loadings")
gq6_parallel <- cfa(
  'f1 =~ gq6_1 + gq6_2 + gq6_3 + gq6_4 + gq6_5',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = c("loadings", "residuals", "lv.variances"))
anova(gq6_congeneric, gq6_tauequivalent, gq6_parallel)
t(sapply(
  list(gq6_congeneric, gq6_tauequivalent, gq6_parallel),
  function(m) fitMeasures(m)[c("chisq", "df", "cfi", "srmr")]
))

## GAC
gac_congeneric <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3',
  data = yg, group = "agegroup", meanstructure = TRUE)
gac_tauequivalent <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = "loadings")
gac_parallel <- cfa(
  'f1 =~ gac_1 + gac_2 + gac_3',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = c("loadings", "residuals", "lv.variances"))
anova(gac_congeneric, gac_tauequivalent, gac_parallel)
t(sapply(
  list(gac_congeneric, gac_tauequivalent, gac_parallel),
  function(m) fitMeasures(m)[c("chisq", "df", "cfi", "srmr")]
))

## GRAT
grat_congeneric <- cfa(
  'f1 =~ losd_2 + losd_3 + losd_4 + losd_5 + losd_6
   f2 =~ sa_1 + sa_2 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_1 + ao_2 + ao_3 + ao_4',
  data = yg, group = "agegroup", meanstructure = TRUE)
grat_tauequivalent <- cfa(
  'f1 =~ losd_2 + losd_3 + losd_4 + losd_5 + losd_6
   f2 =~ sa_1 + sa_2 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_1 + ao_2 + ao_3 + ao_4',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = "loadings")
grat_parallel <- cfa(
  'f1 =~ losd_2 + losd_3 + losd_4 + losd_5 + losd_6
   f2 =~ sa_1 + sa_2 + sa_3 + sa_4 + sa_5 + sa_6
   f3 =~ ao_1 + ao_2 + ao_3 + ao_4',
  data = yg, group = "agegroup", meanstructure = TRUE,
  group.equal = c("loadings", "residuals", "lv.variances"))
anova(grat_congeneric, grat_tauequivalent, grat_parallel)
t(sapply(
  list(grat_congeneric, grat_tauequivalent, grat_parallel),
  function(m) fitMeasures(m)[c("chisq", "df", "cfi", "srmr")]
))
}
}

\keyword{datasets}
