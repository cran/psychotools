\name{plot.raschmodel}

\alias{plot.raschmodel}
\alias{plot.rsmodel}
\alias{plot.pcmodel}
\alias{plot.nplmodel}
\alias{plot.gpcmodel}

\title{Visualizing IRT Models}

\description{
  Base graphics plotting function for various IRT models.
}

\usage{
\S3method{plot}{raschmodel}(x,
  type = c("profile", "curves", "regions", "information", "piplot"), \dots)
\S3method{plot}{rsmodel}(x,
  type = c("regions", "profile", "curves", "information", "piplot"), \dots)
\S3method{plot}{pcmodel}(x,
  type = c("regions", "profile", "curves", "information", "piplot"), \dots)
\S3method{plot}{nplmodel}(x,
  type = c("regions", "profile", "curves", "information", "piplot"), \dots)
\S3method{plot}{gpcmodel}(x,
  type = c("regions", "profile", "curves", "information", "piplot"), \dots)
}

\arguments{
  \item{x}{a fitted model object of class \code{"raschmodel"},
    \code{"rsmodel"}, \code{"pcmodel"}, \code{"nplmodel"} or \code{"gpcmodel"}.}
  \item{type}{character, specifying the type of plot to create.}
  \item{\dots}{further arguments passed over to the specific plotting
    function.}
}

\seealso{\code{\link{curveplot}}, \code{\link{regionplot}},
  \code{\link{profileplot}}, \code{\link{infoplot}}, \code{\link{piplot}}}

\keyword{hplot}

