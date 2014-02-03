\name{plot.PCModel}
\alias{plot.PCModel}


\title{Visualizing Partial Credit Models}


\description{
  Base graphics plotting function for partial credit models.
}


\usage{
\S3method{plot}{PCModel}(x, type = NULL, names = NULL, ref = NULL,
  main = NULL, ylab = "Latent trait", ylim = NULL, off = 0.1, ...)
}


\arguments{
  \item{x}{an object of class \code{"PCModel"}.}

  \item{type}{character, specifying which type of threshold parameters
    are used to mark the category regions per item in the plot (see
    details below).}

  \item{names}{character vector containing x-axis labels for the
  items. If set to \code{NULL} (the default), these are generated.}

  \item{ref}{a vector of labels or position indices of item parameters
    which should be used as restriction/for normalization. If
    \code{NULL} (the default), all items are used (sum zero restriction).}

  \item{main}{character, specifying the overall title of the plot.}

  \item{ylab}{label for the y axis.}

  \item{ylim}{y axis limits.}

  \item{off}{the distance (in scale units) between two item rectangles.}

  \item{\dots}{further arguments passed to \code{\link{plot}}}
}


\details{

  The effect plots visualization implemented here was already used by
  Van der Linden and Hambleton (1997) in the context of IRT and has been
  called "effect plots" by Fox & Hong (2009). In our implementation,
  these plots show, dependent on the choosen type of threshold
  parameters, different regions for the categories of an item over the
  theta axis. If \code{type} is set to \code{"modus"}, the cutpoints
  correspond to the threshold parameters and the rectangles mark the
  theta regions where a category is most probable. If \code{type} is set
  to \code{"median"}, the cutpoints correspond to the point on the theta
  axis, where the cumulative probability to score in category x or
  higher is 0.5. If set to \code{"mean"}, the cutpoints correspond to
  the point on the theta axis where the expected score is exactly
  between two categories, i.e. 0.5 for a dichotomous item.

  If there are unordered threshold parameters, the location of the
  original threshold parameters are indicated by red dashed lines.

}


\references{
  Fox, J. & Hong, Jangman. (2009). Effect Displays in \code{R}
  for Multinomial and Proportional-Odds Logit Models: Extensions to the
  effects Package. \emph{Journal of Statistical Software}, \bold{32(1)},
  \url{http://www.jstatsoft.org/v32/i01/}

  Van der Linden, W. J. & Hambleton, R. K. (1997). \emph{Handbook of
  Modern Item Response Theory}. Springer, New York.
}


\seealso{\code{\link{PCModel.fit}}}%, \code{\link{plotCCC}}}


\examples{

## Verbal aggression data
data("VerbalAggression", package = "psychotools")

## Partial credit model for the first other-to-blame
## situation: 'A bus fails to stop for me'
pcm <- PCModel.fit(VerbalAggression$resp[, 1:6])

## Effect plot with modus as cutpoint and custom labels.
lab <- paste(rep(c("Curse", "Scold", "Shout"), each = 2),
             rep(c("Want", "Do"), 3 ), sep = "-")
plot(pcm, names = lab)

## Compare the cutpoints (with ylim specified manually)
ylim <- c(-2, 2)
layout(matrix(1:3, ncol = 1))
plot(pcm, type = "mode", main = "Modus as Cutpoint", ylim = ylim) 
plot(pcm, type = "median", main = "Median as Cutpoint", ylim = ylim)
plot(pcm, type = "mean", main = "Mean as Cutpoint", ylim = ylim)

## Partial credit model for full VerbalAggression data set
pcm_va <- PCModel.fit(VerbalAggression$resp)
plot(pcm_va)

}


\keyword{hplot}
