\name{print.paircomp}

\alias{print.paircomp}
\alias{format.paircomp}

\title{Formatting Paired Comparison Data}

\description{
  Fine control for formatting and printing objects of \code{"paircomp"} data.
}

\usage{
  \method{format}{paircomp}(x, sep = ", ", brackets = TRUE,
    abbreviate = NULL, width = getOption("width") - 7, \dots)
  \method{print}{paircomp}(x, quote = FALSE, \dots)
}

\arguments{
  \item{x}{an object of class \code{"paircomp"}.}
  \item{sep}{character. A character for separating comparisons within subjects.}
  \item{brackets}{logical or character. Either a logical (Should brackets be
    wrapped around all comparisons for a single subject?) or a character of
    length two with opening and ending symbol.}
  \item{abbreviate}{logical or integer. Should object labels be abbreviated?
    Alternative an integer with the desired abbreviation length. The default
    is some heuristic based on the length of the labels.}
  \item{width}{integer or logical. Maximal width of the string for a subject.
    If \code{FALSE} no maximal width is set.}
  \item{\dots}{arguments passed to other functions.}
  \item{quote}{logical. Should quotes be printed?}
}

\details{
  The \code{print} method just calls \code{format} (passing on all further
  arguments) and then prints the resulting string.
}

\seealso{\code{\link{paircomp}}}

\examples{
pc2 <- paircomp(rbind(
  c(4,  1,  0),
  c(1,  2, -1),
  c(1, -2, -1),
  c(0,  0,  -3)),
  labels = c("New York", "Rio", "Tokyo"))

print(pc2)
print(pc2, abbreviate = FALSE)
print(pc2, abbreviate = FALSE, width = 10)
}

\keyword{classes}
