\name{paircomp}

\alias{paircomp}

\alias{length.paircomp}
\alias{c.paircomp}
\alias{[.paircomp}
\alias{rep.paircomp}
\alias{xtfrm.paircomp}

\alias{as.character.paircomp}
\alias{as.data.frame.paircomp}
\alias{as.double.paircomp}
\alias{as.integer.paircomp}
\alias{as.matrix.paircomp}

\alias{covariates.paircomp}
\alias{covariates<-.paircomp}
\alias{labels.paircomp}
\alias{labels<-.paircomp}
\alias{names.paircomp}
\alias{names<-.paircomp}
\alias{mscale.paircomp}
\alias{mscale<-.paircomp}

\alias{str.paircomp}
\alias{summary.paircomp}
\alias{is.na.paircomp}


\title{Data Structure for Paired Comparisons}

\description{
  A class for representing data from paired comparison experiments
  along with methods for many generic functions.
}

\usage{
  paircomp(data,
    labels = NULL, mscale = NULL, ordered = FALSE, covariates = NULL)
}

\arguments{
  \item{data}{matrix. A matrix with integer values where the rows correspond
    to subjects and the columns to paired comparisons between objects. See
    below for details.}
  \item{labels}{character. A vector of character labels for the objects.
    By default a suitable number of \code{letters} is used.}
  \item{mscale}{integer. A vector of integers giving the measurement scale.
    See below for details. By default guessed from \code{data}.}
  \item{ordered}{logical. Does \code{data} contain both orderings of each
    comparison?}
  \item{covariates}{data.frame. An optional data.frame with object covariates, i.e.,
    it must have the same number of rows as the length of \code{labels}.
    May be \code{NULL} (default).}
}

\details{
  \code{paircomp} is designed for holding paired comparisons of
  \eqn{k} objects measured for \eqn{n} subjects.
  
  The comparisons should be coded in an integer matrix \code{data}
  with \eqn{n} rows (subjects) and \eqn{k \choose 2} columns
  (unless \code{ordered = TRUE}, see below). The columns must be
  ordered so that objects are sequentially compared with all
  previous objects, i.e.: 1:2, 1:3, 2:3, 1:4, 2:4, 3:4, etc.  
  Each column represents the results of a comparison for two particular
  objects. Positive values signal that the first object was preferred,
  negative values that the second was preferred, zero signals no
  preference. Larger absolute values signal stronger preference.
  
  \code{mscale} provides the underlying measurement scale. It must
  be a symmetric sequence of integers of type \code{(-i):i} where
  \code{i} must be at least \code{1}. However, it may exclude
  \code{0} (i.e., forced choice).
  
  If \code{ordered = TRUE}, the order of comparison matters and
  thus \code{data} is assumed to have twice as many columns. The
  second half of columns then corresponds to the comparisons
  2:1, 3:1, 3:2, 4:1, 4:2, 4:3, etc.
}

\value{
  \code{paircomp} returns an object of class \code{"paircomp"} which is
  a matrix (essentially \code{data}) with all remaining arguments
  of \code{paircomp} as attributes (after being
  checked and potentially suitably coerced or transformed).
}


\seealso{\code{\link{subset.paircomp}}, \code{\link{print.paircomp}}}

\examples{
## a simple paired comparison
pc <- paircomp(rbind(
  c(1,  1,  1), # a > b, a > c, b > c
  c(1,  1, -1), # a > b, a > c, b < c
  c(1, -1, -1), # a > b, a < c, b < c
  c(1,  1,  1)))

## basic methods
pc
str(pc)
summary(pc)
pc[2:3]
c(pc[2], pc[c(1, 4)])

## methods to extract/set attributes
labels(pc)
labels(pc) <- c("ah", "be", "ce")
pc
mscale(pc)
covariates(pc)
covariates(pc) <- data.frame(foo = factor(c(1, 2, 2), labels = c("foo", "bar")))
covariates(pc)
names(pc)
names(pc) <- LETTERS[1:4]
pc

## reorder() and subset() both select a subset of
## objects and/or reorders the objects
reorder(pc, c("ce", "ah"))


## include paircomp object in a data.frame
## (i.e., with subject covariates)
dat <- data.frame(
  x = rnorm(4),
  y = factor(c(1, 2, 1, 1), labels = c("hansi", "beppi")))
dat$pc <- pc
dat


## formatting with long(er) labels and extended scale
pc2 <- paircomp(rbind(
  c(4,  1,  0),
  c(1,  2, -1),
  c(1, -2, -1),
  c(0,  0,  -3)),
  labels = c("Nordrhein-Westfalen", "Schleswig-Holstein", "Baden-Wuerttemberg"))
## default: abbreviate
print(pc2)
print(pc2, abbreviate = FALSE)
print(pc2, abbreviate = FALSE, width = FALSE)


## paired comparisons with object covariates
pc3 <- paircomp(rbind(
  c(2,  1,  0),
  c(1,  1, -1),
  c(1, -2, -1),
  c(0,  0,  0)),
  labels = c("New York", "Rio", "Tokyo"),
  covariates = data.frame(hemisphere = factor(c(1, 2, 1), labels = c("North", "South"))))
covariates(pc3)
}

\keyword{classes}
