\name{itemresp}
\alias{itemresp}

\alias{is.na.itemresp}
\alias{labels.itemresp}
\alias{labels<-.itemresp}
\alias{length.itemresp}
\alias{levels.itemresp}
\alias{mscale.itemresp}
\alias{mscale<-.itemresp}
\alias{names.itemresp}
\alias{names<-.itemresp}
\alias{rep.itemresp}
\alias{str.itemresp}
\alias{xtfrm.itemresp}


\title{Data Structure for Item Response Data}

\description{
  A class for representing data from questionnaires
  along with methods for many generic functions.
}

\usage{
  itemresp(data, mscale = NULL, labels = NULL, names = NULL)
}

\arguments{
  \item{data}{matrix or data frame. A matrix or data frame with integer
    values or factors where the rows correspond to subjects and the
    columns to items. See below for details.}
  \item{mscale}{integer or character. A list of vectors (either integer
    or character) giving the measurement scale.
    See below for details. By default guessed from \code{data}.}
  \item{labels}{character. A vector of character labels for the items.
    By default, the column names of \code{data} are used or, if these are
    not available, the string \code{"item"} along with numbers 1, 2, \dots
    is used.}
  \item{names}{character. A vector of names (or IDs) for the subjects. By default,
    no subject names are used.}
}

\details{
  \code{itemresp} is designed for item response data of
  \eqn{n} subjects for \eqn{k} items.
  
  The item responses should be coded in a matrix \code{data}
  with \eqn{n} rows (subjects) and \eqn{k} columns (items). Alternatively,
  \code{data} can be a data frame with \eqn{n} rows (subjects) and
  \eqn{k} variables (items), which can be either factors or integer
  valued vectors.
  
  \code{mscale} provides the underlying measurement scale either as
  integer or character vector(s). If all items are measured on the same
  scale, \code{mscale} can be a vector. Alternatively, it can be
  provided as a named list of vectors for each item. If the list
  contains one unnamed element, this element will be used as the
  measurement scale for items that have not been named. Integers or
  characters not present in \code{mscale} but in \code{data} will be
  replaced by \code{NA}. All items must be measured with at least 2
  categories. By default, \code{mscale} is set to the full range of
  observed values for all integer items (see example below) and the
  corresponding levels for all factor items in \code{data}.

  Methods to standard generic functions include: \code{str},
  \code{length} (number of subjects), \code{dim} (number of subjects and
  items), \code{is.na} (only \code{TRUE} if all item responses are
  \code{NA} for a subject), \code{print} (see
  \code{\link{print.itemresp}} for details), \code{summary} and
  \code{plot} (see \code{\link{summary.itemresp}} for details),
  subsetting via \code{[} and \code{subset} (see
  \code{\link{subset.itemresp}} for details), \code{is.itemresp} and
  various coercion functions to other classes (see
  \code{\link{as.list.itemresp}} for details).
  
  Extracting/replacing properties is available through: \code{labels} for the item labels,
  \code{mscale} for the measurement scale, \code{names} for subject names/IDs.
}

\value{
  \code{itemresp} returns an object of class \code{"itemresp"} which is
  a matrix (\code{data} transformed to integers 0, 1, \dots) plus an
  attribute \code{"mscale"} as a named list for each item
  (after being checked and potentially suitably coerced or transformed
  to all integer or all character).
}

\seealso{\code{\link{print.itemresp}}, \code{\link{summary.itemresp}},
  \code{\link{as.list.itemresp}}, \code{\link{subset.itemresp}}}

\examples{
## binary responses to three items, coded as matrix
x <- cbind(c(1, 0, 1, 0), c(1, 0, 0, 0), c(0, 1, 1, 1))
## transformed to itemresp object
xi <- itemresp(x)

## printing (see also ?print.itemresp)
print(xi)
print(xi, labels = TRUE)

## subsetting/indexing (see also ?subset.itemresp)
xi[2]
xi[c(TRUE, TRUE, FALSE, FALSE)]
subset(xi, items = 1:2)
dim(xi)
length(xi)

## summary/visualization (see also ?summary.itemresp)
summary(xi)
plot(xi)

## query/set measurement scale labels
## extract mscale (tries to collapse to vector)
mscale(xi)
## extract as list
mscale(xi, simplify = FALSE)
## replacement by list
mscale(xi) <- list(item1 = c("no", "yes"),
  item2 = c("nay", "yae"), item3 = c("-", "+"))
xi
mscale(xi)
## replacement with partially named list plus default
mscale(xi) <- list(item1 = c("n", "y"), 0:1)
mscale(xi)
## replacement by vector (if number of categories constant)
mscale(xi) <- c("-", "+")
mscale(xi, simplify = FALSE)

## query/set item labels and subject names
labels(xi)
labels(xi) <- c("i1", "i2", "i3")
names(xi)
names(xi) <- c("John", "Joan", "Jen", "Jim")
print(xi, labels = TRUE)

## coercion (see also ?as.list.itemresp)
## to integer matrix
as.matrix(xi)
## to data frame with single itemresp column
as.data.frame(xi)
## to list of factors
as.list(xi)
## to data frame with factors
as.list(xi, df = TRUE)


## polytomous responses with missing values and a unequal number of
## categories in a data frame
d <- data.frame(q1 = c(-2, 1, -1, 0, NA, 1, NA),
                q2 = c(3, 5, 2, 5, NA, 2, 3),
                q3 = factor(c(1, 2, 1, 2, NA, 3, 2), levels = 1:3,
                            labels = c("disagree", "neutral", "agree")))
di <- itemresp(d)
di

## auto-completion of mscale: full range (-2, ..., 2) for q1, starting
## from smallest observed (negative) value (-2) to the same (positive)
## value (2), full (positive) range for q2, starting from smallest
## observed value (2) to largest observed value (5), missing category of
## 4 is detected, for q3 given factor levels are used
mscale(di)

## set mscale for q2 and add category 1, q1 and q3 are auto-completed:
di <- itemresp(d, mscale = list(q2 = 1:5))

## is.na.itemresp - only true for observation 5 (all missing)
is.na(di)

## illustration for larger data set
data("VerbalAggression", package = "psychotools")
r <- itemresp(VerbalAggression$resp[, 1:12])
str(r)
head(r)
summary(r)
plot(r)


## COPY BACK TO MSACLE.RD WHEN ITEMRESP IS RELEASED.
## similar for "itemresp" data
ir <- itemresp(cbind(
  c(-1, 0, 1, 1, 0),
  c(0, 1, 2, 1, 2),
  c(1, 2, 1, 1, 3)))
ir

## extract
mscale(ir)

## replace (single scale for all items)
mscale(ir) <- 1:3
ir
}

\keyword{classes}
