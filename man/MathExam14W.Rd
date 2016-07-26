\name{MathExam14W}
\alias{MathExam14W}
\encoding{latin1}

\title{Mathematics 101 Exam Results}

\description{
  Responses of 729 students to 13 items in a written
  exam of introductory mathematics along with several covariates. 
}

\usage{data("MathExam14W")}

\format{
A data frame containing 729 observations on 9 variables.
\describe{
  \item{solved}{Item response matrix (of class \code{\link{itemresp}}) with
    values 1/0 coding solved correctly/other.}
  \item{credits}{Item response matrix (of class \code{\link{itemresp}}) with
    values 2/1/0 coding solved correctly/incorrectly/not attempted.}
  \item{nsolved}{Integer. The number of items solved correctly.}
  \item{tests}{Integer. The number of online test exercises solved
    correctly prior to the written exam.}
  \item{gender}{Factor indicating gender.}
  \item{study}{Factor indicating two different types of business/economics
    degrees. Either the 3-year bachelor program (571) or the 4-year
    diploma program (155).}
  \item{semester}{Integer. The number of semesters enrolled in the given
    university program.}
  \item{attempt}{Factor. The number of times the course/exam has
    been attempted (including the current attempt).}
  \item{group}{Factor indicating whether the students were in the
    first or second batch (with somewhat different items) in the exam.}
}
}

\details{
The data provides individual end-term exam results from a Mathematics 101
course for first-year business and economics students at Universit\enc{ä}{ae}t Innsbruck.
The format of the course comprised biweekly online tests (26 numeric exercises,
conducted in OpenOLAT) and a written exam at the end of the semester (13
single-choice exercises with five answer alternatives). The course covers
basics of analysis, linear algebra, financial mathematics, and probability
calculus (where the latter is not assessed in this exam).

In this exam, 729 students participated (out of 941 registered in the course).
To avoid cheating, all students received items with essentially the same questions
but different numbers (using the exams infrastructure of Zeileis et al. 2014).
Also, due to the large number of students two groups of students had to be formed
which received partially different items. The items which differed (namely 1, 5,
6, 7, 8, 9, 11, 12) varied in the setup/story, but not in the 
mathematical skills needed to solve the exercises. Prior to the exam,
the students could select themselves either into the first group (early in the morning)
or the second group (starting immediately after the end of the first group).

Correctly solved items yield 100 percent of the associated points. Items
without correct solution can either be unanswered (0 percent) or receive an
incorrect answer (minus 25 percent) to discourage random guessing. In the examples
below, the items are mostly only considered as binary. Typically, students with
8 out of 13 correct answers passed the course.
}

\source{
Department of Statistics, Universit\enc{ä}{ae}t Innsbruck
}

\references{
Zeileis A, Umlauf N, Leisch F (2014).
  Flexible Generation of E-Learning Exams in R: Moodle Quizzes, OLAT Assessments, and Beyond.
  \emph{Journal of Statistical Software}, \bold{58}(1), 1--36.
  doi:10.18637/jss.v058.i01
}

\seealso{\code{\link{itemresp}}, \code{\link{raschmodel}}, \code{\link{pcmodel}}, \code{\link{anchortest}}}

\examples{
## load data and exclude extreme scorers
data("MathExam14W", package = "psychotools")
me <- subset(MathExam14W, nsolved > 0 & nsolved < 13)

## item response data:
## solved (correct/other) or credits (correct/incorrect/not attempted)
par(mfrow = c(1, 2))
plot(me$solved)
plot(me$credits)


## Rasch and partial credit model
ram <- raschmodel(me$solved)
pcm <- pcmodel(me$credits)

## various types of graphics displays
plot(ram, type = "profile")
plot(pcm, type = "profile", add = TRUE, col = "blue")
plot(ram, type = "piplot")
plot(pcm, type = "piplot")
plot(ram, type = "region")
plot(ram, type = "curves")
plot(pcm, type = "curves")


\donttest{
## test for differential item function with automatic anchoring
## passing vs. not passing students
at1 <- anchortest(solved ~ factor(nsolved <= 7), data = me,
  adjust = "single-step")
at1
plot(at1$final_tests)
## -> "good" students discriminate somewhat more
## (quad/payflow/lagrange are slightly more difficult)

## group 1 vs. group 2
at2 <- anchortest(solved ~ group, data = me, adjust = "single-step")
at2
plot(at2$final_tests)
## -> quad/payflow/planning easier for group 1
## -> interest slightly easier for group 2

## bring out differences between groups 1 and 2
## by (anchored) item difficulty profiles
ram1 <- raschmodel(subset(me, group == "1")$solved)
ram2 <- raschmodel(subset(me, group == "2")$solved)
plot(ram1, parg = list(ref = at2$anchor_items), ylim = c(-2, 3))
plot(ram2, parg = list(ref = at2$anchor_items), add = TRUE, col = "blue")
legend("topleft", c("Group 1", "Group 2"), pch = 21,
  pt.bg = c("lightgray", "blue"), bty = "n")
}
}

\keyword{datasets}
