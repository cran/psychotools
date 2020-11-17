\name{plot.paircomp}
\alias{plot.paircomp}

\title{Plotting Paired Comparison Data}

\description{
  Plotting the frequency table from \code{"paircomp"} data.
}

\usage{
  \method{plot}{paircomp}(x, off = 0.05,
  xlab = "Proportion of comparisons", ylab = "", tol.xlab = 0.05,
  abbreviate = TRUE, hue = NULL, chroma = 40, luminance = 80,
  xlim = c(0, 1), ylim = NULL, xaxs = "i", yaxs = "i", \dots)
}

\arguments{
  \item{x}{an object of class \code{"paircomp"}.}
  \item{off}{numeric. Offset between segments on the y-axis.}
  \item{xlab, ylab}{character. Axis labels.}
  \item{tol.xlab}{numeric. convenience tolerance parameter for x-axis annotation.
    If the distance between two labels drops under this threshold, they
    are plotted equidistantly.}
  \item{abbreviate}{logical or integer. Should object labels be abbreviated?
    Alternative an integer with the desired abbreviation length. The default
    is some heuristic based on the length of the labels.}
  \item{hue}{numeric. A vector of hues in [0, 360], recycled to the
    number of objects compared in \code{x}. A sequential palette
    is computed for each hue, see below.}
  \item{chroma}{numeric. Maximum chroma in the palette.}
  \item{luminance}{numeric. Minimum (and maximum) luminance in the
    palette. If omitted, the maximum is set to 95.}
  \item{xlim, ylim, xaxs, yaxs, \dots}{graphical arguments passed to
    \code{\link{plot}}.}
}

\details{
  The \code{plot} method creates a frequency table (using \code{summary})
  and visualizes this using a sort of spine plot with HCL-based
  diverging palettes. See Zeileis, Hornik, Murrell (2009) for the
  underlying ideas.
}

\references{
  Zeileis A, Hornik K, Murrell P (2009),
  Escaping RGBland: Selecting Colors for Statistical Graphics.
  \emph{Computational Statistics & Data Analysis}, \bold{53}, 3259-3270.
  \doi{10.1016/j.csda.2008.11.033}
}

\seealso{\code{\link{paircomp}}}

\examples{
data("GermanParties2009", package = "psychotools")
par(mar = c(5, 6, 3, 6))
plot(GermanParties2009$preference, abbreviate = FALSE)
}

\keyword{classes}
