\name{subset.paircomp}

\alias{subset.paircomp}
\alias{reorder.paircomp}

\title{Subsetting/Reordering Paired Comparison Data}

\description{
  Selection of subsets of objects to be compared and/or reordering of  
  objects in \code{"paircomp"} data.
}

\usage{
  \method{reorder}{paircomp}(x, labels, \dots)
  \method{subset}{paircomp}(x, subset, select, \dots)

}

\arguments{
  \item{x}{an object of class \code{"paircomp"}.}
  \item{labels, select}{character or integer. Either a vector of
    (at least two) elements of \code{labels(x)} or an integer
    with their position. Partial string matching is enabled.}
  \item{subset}{currently not implemented. (Should be a specification
    of subsets of subjects.)}
  \item{\dots}{currently not used.}
}

\details{
  The \code{subset} method currently just calls the \code{reorder} method.
}

\seealso{\code{\link{paircomp}}}

\examples{
pc <- paircomp(rbind(
  c(1,  1,  1), # a > b, a > c, b > c
  c(1,  1, -1), # a > b, a > c, b < c
  c(1, -1, -1), # a > b, a < c, b < c
  c(1,  1,  1)))
reorder(pc, c("c", "a"))
}

\keyword{classes}
