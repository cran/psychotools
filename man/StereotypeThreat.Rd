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
## Step 2: Fix loadings across groups
m2  <- c(
'ability   =~ label(rep("load_n", 4)) * numerical  + 1 * abstract + label(rep("load_v", 4)) * verbal',
'ability   ~  0 * 1',
'ability   ~~ label(c("lvar:maj_c", "lvar:min_c", "lvar:maj_t", "lvar:min_t")) * ability',

'numerical ~  label(c("mean_n:maj_c", "mean_n:min_c", "mean_n:maj_t", "mean_n:min_t")) * 1',
'abstract  ~  label(c("mean_a:maj_c", "mean_a:min_c", "mean_a:maj_t", "mean_a:min_t")) * 1',
'verbal    ~  label(c("mean_v:maj_c", "mean_v:min_c", "mean_v:maj_t", "mean_v:min_t")) * 1',

'numerical ~~ label(c("var_n:maj_c", "var_n:min_c", "var_n:maj_t", "var_n:min_t")) * numerical',
'abstract  ~~ label(c("var_a:maj_c", "var_a:min_c", "var_a:maj_t", "var_a:min_t")) * abstract',
'verbal    ~~ label(c("var_v:maj_c", "var_v:min_c", "var_v:maj_t", "var_v:min_t")) * verbal')

## Step 2a: Free numerical loading in group min_t (minority.threat)
m2a <- m2
m2a[1] <- 'ability   =~ label(c(rep("load_n", 3), "load_n:min_t")) * numerical  + 1 * abstract + label(rep("load_v", 4)) * verbal'

## Step 3: Fix variances across groups
m3 <- m2a
m3[7:9] <- c(
'numerical ~~ label(rep("var_n", 4)) * numerical',
'abstract  ~~ label(rep("var_a", 4)) * abstract',
'verbal    ~~ label(rep("var_v", 4)) * verbal')

## Step 3a: Free numerical variance in group min_t
m3a <- m3
m3a[7] <- 'numerical ~~ label(c(rep("var_n", 3), "var_n:min_t")) * numerical'

## Step 4: Fix latent variances within conditions
m4 <- m3a
m4[3] <- 'ability   ~~ label(c("lvar:maj", "lvar:min", "lvar:maj", "lvar:min")) * ability'

## Step 5: Fix certain means, free others
m5 <- m4
m5[c(2, 4:6)] <- c(
'ability   ~  label(c(NA, "lmean:min_c", "lmean:maj_t", "lmean:min_t")) * 1 + c(0, NA, NA, NA) * 1',
'numerical ~  label(c(rep("mean_n", 3), "mean_n:min_t")) * 1',
'abstract  ~  label(rep("mean_a", 4)) * 1',
'verbal    ~  label(rep("mean_v", 4)) * 1')

## Step 5a: Free ability mean in group maj_c
m5a <- m5
m5a[5] <- 'abstract  ~  label(c("mean_a:maj_c", rep("mean_a", 3))) * 1'

## Step 5b: Free ability mean in group min_c
m5b <- m5a
m5b[5] <- 'abstract  ~  label(c("mean_a:maj_c", "mean_a:min_c", rep("mean_a", 2))) * 1'

## Step 6: Different latent mean structure
m6 <- m5b
m6[2] <- 'ability   ~  label(c("lmean:maj", NA, "lmean:maj", NA)) * 1 + c(NA, 0, NA, 0) * 1'

## Fit all models
m <- list(m2, m2a, m3, m3a, m4, m5, m5a, m5b, m6)
names(m) <- c("2", "2a", "3", "3a", "4", "5", "5a", "5b", "6")
fm <- lapply(m, function(f)
  lavaan(f, data = StereotypeThreat, meanstructure = TRUE, group = "group", likelihood = "wishart"))

## Extract measures of fit
tab <- t(sapply(fm, fitMeasures, c("chisq", "df", "pvalue", "rmsea", "cfi")))
tab <- rbind("1" = c(0, 0, 1, 0, 1), tab)
tab <- cbind(tab,
  delta_chisq = c(NA, abs(diff(tab[, "chisq"]))),
  delta_df = c(NA, diff(tab[, "df"])))
tab <- cbind(tab, "pvalue2" = pchisq(tab[, "delta_chisq"], abs(tab[, "delta_df"]), lower.tail = FALSE))
tab <- tab[, c(2, 1, 3, 7, 6, 8, 4, 5)]
round(tab, digits = 3)
}
}

\keyword{datasets}
