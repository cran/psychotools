\name{StereotypeThreat}
\alias{StereotypeThreat}

\title{Stereotype Threat in Dutch Differential Aptitude Test}

\description{
  Cross-section data from Differential Aptitude Test (DAT) among Dutch highschool
  students, along with experimental conditions pertaining to stereotype threat.
}

\usage{data("StereotypeThreat")}

\format{
A data frame containing 295 observations on 11 variables.
\describe{
  \item{condition}{Factor indicating experimental condition: \code{"control"} or
    stereotype \code{"threat"}, for details see below.}
  \item{ethnicity}{Factor coding ethnicity: Dutch \code{"majority"} or \code{"minority"}.}
  \item{numerical}{Number of items solved in numerical ability subtest
    (out of 14 complicated mathematical items).}
  \item{abstract}{Number of items solved in abstract reasoning subtest
    (out of 18 items with a logical sequence of diagrams).}
  \item{verbal}{Number of items solved in verbal reasoning subtest
    (out of 16 verbal analogy items).}
  \item{gender}{Factor indicating gender.}
  \item{age}{Age in years.}
  \item{vintelligence}{Numerical coding of the value of one's own intelligence.
    Answer to: How important is your intelligence for you?
    Range is from very important (\code{5}) to unimportant (\code{1}).}
  \item{vgrades}{Numerical coding of the value of getting good grades.  
    Answer to: How much do you value getting good school grades?
    Range is from a lot of value (\code{5}) to not so much value (\code{1}).}
  \item{vprejudice}{Numerical coding of the answer to: Do you think that people of your group are prejudiced against?
    Range is from certainly (\code{5}) to not at all (\code{1}).}
  \item{gpa}{Numerical grade point average on 10-point scale (with 10 being the best grade).
    It has 57 missing values as some schools were either unwilling to share the
    data or did not provide it timely enough.}
}
}

\details{
The data are taken from Study 1 of Wicherts et al. (2005) and have been used to study
stereotype threat on intelligence test performance among Dutch highschool students.

On average, Dutch minority students attain lower educational levels compared to Dutch
majority students and studies have shown that minority students are often viewed as
less smart/educated. Conversely, minorities often feel discriminated against in scholastic
domains.

Wicherts et al. (2005) administered an intelligence test consisting of three
subtests (for numerical ability, abstract reasoning, and verbal reasoning) and varied
the amount of stereotype threat related to ethnic minorities by changing the presentation
of the test. In the \code{"threat"} condition, the questions were declared to be part
of an intelligence test and also an ethnicity questionnaire was conducted prior to the DAT.
In the \code{"control"} condition, intelligence was not mentioned and no ethnicity
questionnaire was conducted.

The variables \code{numerical}, \code{abstract}, and \code{verbal} can be used to assess
ability/intelligence. And the \code{vintelligence}, \code{vgrades}, \code{vprejudice}, and
\code{gpa} variables capture identification with the scholastic domain.

See Wicherts et al. (2005) for details.
}

\source{
Provided by Jelte M. Wicherts.
}

\references{
Wicherts, J.M., Conor, V.D., Hessen, D.J. (2005).
  Stereotype Threat and Group Differences in Test Performance: A Question of Measurement Invariance.
  \emph{Journal of Personality and Social Psychology}, \bold{89}(5), 696-716.
}

\examples{
## Data: Load and include/order wrt group variable
data("StereotypeThreat", package = "psychotools")
StereotypeThreat <- transform(StereotypeThreat, group = interaction(ethnicity, condition))
StereotypeThreat <- StereotypeThreat[order(StereotypeThreat$group),]

## Exploratory analysis (Table 2, p. 703)
tab2 <- with(StereotypeThreat, rbind(
   "#"         = tapply(numerical, group, length),
   "Numerical" = tapply(numerical, group, mean),
   "         " = tapply(numerical, group, sd),
   "Abstract " = tapply(abstract,  group, mean),
   "         " = tapply(abstract,  group, sd),
   "Verbal   " = tapply(verbal,    group, mean),
   "         " = tapply(verbal,    group, sd)))
round(tab2, digits = 2)

## Corresponding boxplots
plot(numerical ~ group, data = StereotypeThreat)
plot(abstract  ~ group, data = StereotypeThreat)
plot(verbal    ~ group, data = StereotypeThreat)

## MANOVA (p. 703)
m <- lm(cbind(numerical, abstract, verbal) ~ ethnicity * condition, data = StereotypeThreat)
anova(m, update(m, . ~ . - ethnicity:condition))
## corresponding univariate results
printCoefmat(t(sapply(summary(m),
  function(x) x$coefficients["ethnicityminority:conditionthreat", ])))

## MGCFA (Table 3, p. 704)
## can be replicated using package lavaan
\dontrun{
## convenience function for multi-group CFA on this data
mgcfa <- function(model, ...) cfa(model, data = StereotypeThreat,
  group = "group", likelihood = "wishart", start = "simple", ...)
		 
## list of all 9 models
m <- vector("list", length = 9)
names(m) <- c("m2", "m2a", "m3", "m3a", "m4", "m5", "m5a", "m5b", "m6")

## Step 2: Fix loadings across groups
f <- 'ability =~ abstract + verbal + numerical'
m$m2 <- mgcfa(f, group.equal = "loadings")

## Step 2a: Free numerical loading in group 4 (minority.threat)
f <- 'ability =~ abstract + verbal + c(l1, l1, l1, l4) * numerical'
m$m2a <- mgcfa(f, group.equal = "loadings")

## Step 3: Fix variances across groups
m$m3 <- mgcfa(f, group.equal = c("loadings", "residuals"))

## Step 3a: Free numerical variance in group 4
f <- c(f, 'numerical ~~ c(e1, e1, e1, e4) * numerical')
m$m3a <- mgcfa(f, group.equal = c("loadings", "residuals"))

## Step 4: Fix latent variances within conditions
f <- c(f, 'ability ~~ c(vmaj, vmin, vmaj, vmin) * ability')
m$m4 <- mgcfa(f, group.equal = c("loadings", "residuals"))

## Step 5: Fix certain means, free others
f <- c(f, 'numerical ~ c(na1, na1, na1, na4) * 1')
m$m5 <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))

## Step 5a: Free ability mean in group majority.control
f <- c(f, 'abstract ~ c(ar1, ar2, ar2, ar2) * 1')
m$m5a <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))

## Step 5b: Free also ability mean in group minority.control
f <- c(f[1:4], 'abstract ~ c(ar1, ar2, ar3, ar3) * 1')
m$m5b <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))

## Step 6: Different latent mean structure
f <- c(f, 'ability ~  c(maj, min1, maj, min2) * 1 + c(NA, 0, NA, 0) * 1')
m$m6 <- mgcfa(f, group.equal = c("loadings", "residuals", "intercepts"))

## Extract measures of fit
tab <- t(sapply(m, fitMeasures, c("chisq", "df", "pvalue", "rmsea", "cfi")))
tab <- rbind("1" = c(0, 0, 1, 0, 1), tab)
tab <- cbind(tab,
  delta_chisq = c(NA, abs(diff(tab[, "chisq"]))),
  delta_df = c(NA, diff(tab[, "df"])))
tab <- cbind(tab, "pvalue2" = pchisq(tab[, "delta_chisq"],
  abs(tab[, "delta_df"]), lower.tail = FALSE))
tab <- tab[, c(2, 1, 3, 7, 6, 8, 4, 5)]
round(tab, digits = 3)
}
}

\keyword{datasets}
