\name{summary.itemresp}

\alias{summary.itemresp}
\alias{plot.itemresp}

\title{Summarizing and Visualizing Item Response Data}

\description{
  Summarizing and visualizing \code{"itemresp"} data objects.
}

\usage{
  \method{summary}{itemresp}(object, items = NULL, abbreviate = FALSE,
    mscale = TRUE, simplify = TRUE, sep = " ", \dots)
  \method{plot}{itemresp}(x, xlab = "", ylab = "", items = NULL,
    abbreviate = FALSE, mscale = TRUE, sep = "\n", off = 2, axes = TRUE,
    names = TRUE, srt = 45, adj = c(1.1, 1.1), \dots)
}

\arguments{
  \item{object, x}{an object of class \code{"itemresp"}.}
  \item{items}{character or integer for subsetting the items to be
    summarized/visualized. By default, all items are used.}
  \item{abbreviate}{logical or integer. Should scale labels be abbreviated?
    Alternatively, an integer with the desired abbreviation length. The default
    is some heuristic based on the length of the labels.}
  \item{mscale}{logical. Should mscale values be used for printing/plotting?
    If \code{FALSE}, integers 0, 1, \dots are used.}
  \item{simplify}{logical. Should the summary table be collapsed into
    a matrix or returned as a list?}
  \item{sep}{character. A character for separating item labels from
    their corresponding scale labels (if any).}   
  \item{xlab, ylab, off, axes, \dots}{arguments passed to \code{\link[graphics]{spineplot}}.}
  \item{names}{logical or character. If \code{TRUE}, the names of 
    the items are displayed on the x-axis. If \code{FALSE}, numbers of
    items are shown. Alternatively a character vector of the same
    length as the number of items can be supplied.}
  \item{srt, adj}{numeric. Angle (\code{srt}) and adjustment (\code{adj})
    in case names (rather than numbers) are used as x-axis labels.
    These are passed to \code{\link[graphics]{text}}.}
}

\details{
  The \code{plot} method essentially just calls \code{summary} (passing on most further
  arguments) and then visualizes the result as a \code{spineplot}.
}

\seealso{\code{\link{itemresp}}, \code{\link[graphics]{spineplot}}}

\examples{
## summary/visualization for verbal aggression data
data("VerbalAggression", package = "psychotools")
r <- itemresp(VerbalAggression$resp[, 1:6])
mscale(r) <- c("no", "perhaps", "yes")
summary(r)
plot(r)

## modify formatting of mscale
summary(r, abbreviate = 1)
summary(r, mscale = FALSE)

## illustration for varying mscale across items
## merge with additional random binary response
b <- itemresp(rep(c(-1, 1), length.out = length(r)),
  mscale = c(-1, 1), labels = "Dummy")
rb <- merge(r[, 1:2], b)
head(rb, 2)
## summary has NAs for non-existent response categories
summary(rb)
summary(rb, mscale = FALSE)
plot(rb, srt = 25)
plot(rb, mscale = FALSE)
}

\keyword{classes}
