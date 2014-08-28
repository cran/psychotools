\name{subset.itemresp}

\alias{subset.itemresp}
\alias{merge.itemresp}
\alias{c.itemresp}
\alias{[.itemresp}

\title{Subsetting Item Response Data}

\description{
  Subsetting and combining \code{"itemresp"} data objects.
}

\usage{
  \method{subset}{itemresp}(x, items = NULL, subjects = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{"itemresp"}.}
  \item{items}{character, integer, or logical for subsetting the items.}
  \item{subjects}{character, integer, or logical for subsetting the subjects.}
   \item{\dots}{currently not used.}
}

\details{
  The \code{subset} method selects subsets of items and/or subjects in
  item response data. Alternatively, the \code{[} method can be used
  with the row index corresponding to subjects and the column index
  corresponding to items. 
  
  The \code{c} method can be used to combine item response data from
  different subjects for the same items

  The \code{merge} method can be used to combine item response data
  from the same subjects for different items.
}

\seealso{\code{\link{itemresp}}}

\examples{
## binary responses to three items, coded as matrix
x <- cbind(c(1, 0, 1, 0), c(1, 0, 0, 0), c(0, 1, 1, 1))
xi <- itemresp(x)

## subsetting/indexing
xi[2]
xi[-(3:4)]
xi[c(TRUE, TRUE, FALSE, FALSE)]
subset(xi, items = 1:2) # or xi[, 1:2]
subset(xi, items = -2, subjects = 2:3)

## combine two itemresp vectors for different subjects but the same items
xi12 <- xi[1:2]
xi34 <- xi[3:4]
c(xi12, xi34)

## combine two itemresp vectors for the same subjects but different items
## polytomous responses in a data frame
d <- data.frame(q1 = c(-2, 1, -1, 0), q2 = factor(c(1, 3, 1, 3),
  levels = 1:3, labels = c("disagree", "neutral", "agree")))
di <-itemresp(d)
merge(xi, di)

## if subjects have names/IDs, these are used for merging
names(xi) <- c("John", "Joan", "Jen", "Jim")
names(di) <-         c("Joan", "Jen", "Jim", "Jo")
merge(xi, di)
merge(xi, di, all = TRUE)
}

\keyword{classes}
