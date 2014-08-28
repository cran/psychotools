\name{print.itemresp}

\alias{print.itemresp}
\alias{format.itemresp}

\title{Formatting Item Response Data}

\description{
  Fine control for formatting and printing \code{"itemresp"} data objects.
}

\usage{
  \method{format}{itemresp}(x, sep = c(",", ":"), brackets = TRUE,
    abbreviate = NULL, mscale = TRUE, labels = FALSE,
    width = getOption("width") - 7L, \dots)
  \method{print}{itemresp}(x, quote = FALSE, \dots)
}

\arguments{
  \item{x}{an object of class \code{"itemresp"}.}
  \item{sep}{character. A character of length 2 (otherwise expanded/reduced)
    for separating responses and items, respectively.}
  \item{brackets}{logical or character. Either a logical (Should brackets be
    wrapped around all responses for a single subject?) or a character of
    length 2 with opening and ending symbol.}
  \item{abbreviate}{logical or integer. Should scale labels be abbreviated?
    Alternatively, an integer with the desired abbreviation length. The default
    is some heuristic based on the length of the labels.}
  \item{mscale}{logical. Should mscale values be used for printing?
    If \code{FALSE}, integers 0, 1, \dots are used.}
  \item{labels}{logical. Should item labels be displayed?}
  \item{width}{integer or logical. Maximal width of the string for a subject.
    If \code{FALSE} no maximal width is set.}
  \item{\dots}{arguments passed to other functions.}
  \item{quote}{logical. Should quotes be printed?}
}

\details{
  The \code{print} method just calls \code{format} (passing on all further
  arguments) and then prints the resulting string.
}

\seealso{\code{\link{itemresp}}}

\examples{
## item responses from binary matrix
x <- cbind(c(1, 0, 1, 0), c(1, 0, 0, 0), c(0, 1, 1, 1))
xi <- itemresp(x)
## change mscale
mscale(xi) <- c("-", "+")
xi

## flexible formatting
## no/other brackets
print(xi, brackets = FALSE)
print(xi, brackets = c(">>", "<<"))

## include item labels (with different separators)
print(xi, labels = TRUE)
print(xi, labels = TRUE, sep = c(" | ", ": "))

## handling longer mscale categories
mscale(xi) <- c("disagree", "agree")
print(xi)
print(xi, mscale = FALSE)
print(xi, abbreviate = FALSE)
print(xi, abbreviate = FALSE, width = 23)
print(xi, abbreviate = 2)
}

\keyword{classes}
