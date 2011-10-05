\name{SoundQuality}
\alias{SoundQuality}

\title{Quality of Multichannel Reproduced Sound}

\description{
  Paired comparison judgments of 40 selected listeners with respect to
  eight audio reproduction modes and for types of music. 
}

\usage{data("SoundQuality")}

\format{
A data frame containing 783 observations on 6 variables.
\describe{
  \item{id}{Factor. Listener ID.}
  \item{time}{Factor. Listening experiment before or after elicitation and
    scaling of more specific auditory attributes.}
  \item{progmat}{Factor. The Program material: Beethoven, Rachmaninov, Steely
    Dan, Sting.}
  \item{repet}{The repetition within each time point.}
  \item{session}{The experimental session coding the presentation order of
    the program material.}
  \item{preference}{Paired comparison of class \code{\link{paircomp}}.
    Preferences for all 28 paired comparisons from 8 audio reproduction modes:
    Mono, Phantom Mono, Stereo, Wide-Angle Stereo, 4-channel Matrix,
    5-channel Upmix 1, 5-channel Upmix 2, and 5-channel Original.}
}
}

\details{
  The data were collected within a series of experiments conducted at the
  Sound Quality Research Unit (SQRU), Department of Acoustics, Aalborg
  University, Denmark, between September 2004 and March 2005.

  The results of scaling listener preference and spatial and timbral auditory
  attributes are reported in Choisel and Wickelmaier (2007).

  Details about the loudspeaker setup and calibration are given in Choisel and
  Wickelmaier (2006).

  The attribute elicitation procedure is described in Wickelmaier and
  Ellermeier (2007) and in Choisel and Wickelmaier (2006).

  The selection of listeners for the experiments is described in Wickelmaier
  and Choisel (2005).
}

\references{
  Choisel, S. & Wickelmaier, F. (2006). Extraction of auditory features and
  elicitation of attributes for the assessment of multichannel reproduced
  sound. Journal of the Audio Engineering Society, 54(9), 815-826.

  Choisel, S. & Wickelmaier, F. (2007). Evaluation of multichannel
  reproduced sound: Scaling auditory attributes underlying listener
  preference. Journal of the Acoustical Society of America, 121(1),
  388-400. 

  Wickelmaier, F. & Choisel, S. (2005). Selecting participants for listening
  tests of multichannel reproduced sound. Presented at the AES 118th
  Convention, May 28-31, Barcelona, Spain, convention paper 6483.

  Wickelmaier, F. & Ellermeier, W. (2007). Deriving auditory features from
  triadic comparisons. Perception & Psychophysics, 69(2), 287-297.
}

\seealso{\code{\link{paircomp}}}

\examples{
data("SoundQuality", package = "psychotools")
summary(SoundQuality$preference)
ftable(xtabs(~ time + repet + progmat, data = SoundQuality))
}

\keyword{datasets}
