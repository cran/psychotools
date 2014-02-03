\name{plot.RSModel}
\alias{plot.RSModel}


\title{Visualizing Rating Scale Models}


\description{
  Base graphics plotting function for rating scale models.
}


\usage{
\S3method{plot}{RSModel}(x, pattern = TRUE, names = NULL, ref = NULL,
  main = NULL, ylab = "Latent trait", ylim = NULL,
  col = c("gray", "gray"), pch = c(21, 22), lty = c(2, 3),
  refline = TRUE, reflinecol = "lightgray", ...)
}


\arguments{
  \item{x}{an object of class \code{"RSModel"}.}

  \item{pattern}{logical. If \code{TRUE}, the pattern of the item difficulties and
    threshold parameters is drawn, otherwise an effect plot similar to
    that of \code{\link{plot.PCModel}} is drawn, where additional
    arguments can be specified via the \code{\dots} argument (see help
    page of \code{\link{plot.PCModel}} for additional informations)}

  \item{names}{character vector containing x-axis labels for the item- and
    treshold parameters. If set to \code{NULL} (the default), these are generated.}

  \item{ref}{logical. If \code{TRUE} (the default), the threshold
    parameters are normed to sum zero, otherwise (\code{FALSE}), the first
    threshold parameter is set to zero.}

  \item{main}{character, specifying the overall title of the plot.}

  \item{ylab}{label for the y axis.}

  \item{ylim}{y axis limits}

  \item{col}{character vector of length two with color names for filling
  the points of the item difficulty and threshold parameters.}

  \item{pch}{numeric vector of length two, specifying the symbols used
    for drawing the item difficulty and threshold parameters (see
    \code{\link[graphics]{points}} for possible values).}
  
  \item{lty}{numeric vector of length two, specifying the line type used
    for the connection lines between the various item difficulty and
    threshold parameters (see \code{\link[graphics]{lines}} for possible values).}
  
  \item{refline}{logical. Should a horizontal line for the reference level
    (mean item difficulty) be drawn?}

  \item{reflinecol}{line color for reference line (if \code{ref}).}
  
  \item{\dots}{further arguments passed to \code{\link[graphics]{plot}} or if
    \code{pattern} is \code{TRUE}, arguments similar to that of \code{\link{plot.PCModel}}.}
}


\seealso{\code{\link{RSModel.fit}}, \code{\link{plot.PCModel}}}%, \code{\link{plotCCC}}}


\examples{

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Rating scale model for the first other-to-blame
## situation: 'A bus fails to stop for me'
rsm <- RSModel.fit(VerbalAggression$resp[, 1:6])

## Basic plot with parameter patterns and custom labels.
lab <- c(paste(rep(c("Curse", "Scold", "Shout"), each = 2),
             rep(c("Want", "Do"), 3 ), sep = "-"), "Tau 1", "Tau 2")
plot(rsm, names = lab)

## Basic plot again, but omitting reference and connection lines
plot(rsm, refline = FALSE, lty = c(0, 0))

## Effects plot, similar to effect plots for PCModel's, with
## 'modus' as cutpoint via additional argument 'type'.
plot(rsm, pattern = FALSE, type = "mode")

}


\keyword{hplot}
