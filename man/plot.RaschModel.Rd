\name{plot.RaschModel}
\alias{plot.RaschModel}

\title{Visualizing Simple Rasch Models}

\description{
  Base graphics plotting function for Rasch models.
}

\usage{
\S3method{plot}{RaschModel}(x, difficulty = TRUE,
  center = TRUE, index = TRUE, names = NULL, abbreviate = FALSE, ref = TRUE,
  col = cbind("lightgray", "black"), refcol = "lightgray", linecol = "black", lty = 2,
  cex = 1, pch = cbind(19, 1), type = NULL, ylim = NULL, xlab = "Items", ylab = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{"RaschModel"}.}
  \item{difficulty}{logical. Should item difficulty (or alternatively: easiness)
    parameters be displayed?}
  \item{center}{logical. Should the item parameters be centered?}
  \item{index}{logical. Should different indexes for different items be used?}
  \item{names}{logical. Should the names for the objects be displayed?}
  \item{abbreviate}{logical or numeric. Should object names be abbreviated?
    If numeric this controls the length of the abbreviation.}
  \item{ref}{logical. Should a horizontal line for the reference level be drawn?}
  \item{col, pch, cex}{graphical appearance of plotting symbols. Can be of the same
    length as number of items. Additionally \code{col} and \code{pch} can
    be matrices with two columns resulting in two symbols being overplotted.}
  \item{refcol}{line color for reference line (if \code{ref}).}
  \item{linecol}{line color.}
  \item{lty}{line type.}
  \item{type}{plot type. Default is \code{"b"} if \code{index} is \code{TRUE}.}
  \item{ylim}{y axis limits.}
  \item{xlab, ylab}{x and y axis labels.}
  \item{\dots}{further arguments passed to \code{\link{plot}}.}
}

\seealso{\code{\link{RaschModel.fit}}}

\keyword{hplot}
