\name{as.list.itemresp}

\alias{is.itemresp}
\alias{as.list.itemresp}
\alias{as.character.itemresp}
\alias{as.data.frame.itemresp}
\alias{as.double.itemresp}
\alias{as.integer.itemresp}
\alias{as.matrix.itemresp}

\title{Coercing Item Response Data}

\description{
  Coercing \code{"itemresp"} data objects to other classes.
}

\usage{
  \method{as.list}{itemresp}(x, items = NULL, mscale = TRUE, df = FALSE, \dots)
}

\arguments{
  \item{x}{an object of class \code{"itemresp"}.}
  \item{items}{character, integer, or logical for subsetting the items.}
  \item{mscale}{logical. Should the measurement scale labels be used
    for creating factor levels? If \code{FALSE}, the values 0, 1, \dots
    are used.}
  \item{df}{logical. Should a data frame of factors be returned?
    If \code{FALSE}, a plain list of factors is returned.}
  \item{\dots}{currently not used.}
}

\details{
  The \code{as.list} method coerces item response data to a list
  (or data frame) of factors with factor levels either taken from
  the \code{mscale(x)} or as the values 0, 1, \dots.
  
  The \code{as.data.frame} method returns a data frame with a single
  column of class \code{"itemresp"}.
  
  Furthermore, \code{as.matrix}, \code{as.integer}, \code{as.double}
  all return a matrix with the item responses coded as values 0, 1, \dots
  
  The \code{as.character} method simply calls \code{\link{format.itemresp}}.

  \code{is.itemresp} can be used to check wether a given object is of
  class \code{"itemresp"}.
}

\seealso{\code{\link{itemresp}}}

\examples{
## item responses from binary matrix
x <- cbind(c(1, 0, 1, 0), c(1, 0, 0, 0), c(0, 1, 1, 1))
xi <- itemresp(x)
## change mscale
mscale(xi) <- c("-", "+")
xi

## coercion to list of factors with levels taken from mscale
as.list(xi)
## same but levels taken as integers 0, 1
as.list(xi, mscale = FALSE)
## only for first two items
as.list(xi, items = 1:2)
## result as data.frame
as.list(xi, df = TRUE)

## data frame with single itemresp column
as.data.frame(xi)

## integer matrix
as.matrix(xi)

## character vector
as.character(xi)

## check class of xi
is.itemresp(xi)
}

\keyword{classes}
