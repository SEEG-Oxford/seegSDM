\name{occurrence2SPDF}
\alias{occurrence2SPDF}

\title{
Convert occurrence dataframe to a SpatialPointsDataFrame
}
\description{
A simple helper function to convert a datframe containing occurrence data into a \code{SpatialPointsDataFrame} for use with other \code{seegSDM} functions.
}
\usage{
occurrence2SPDF(occurrence, crs=wgs84(TRUE))
}

\arguments{
  \item{occurrence}{
A dataframe containing \emph{at least} the columns \code{Latitude} and \code{Longitude} coordinates.
}
 \item{crs}{
The CRS for the occurrence data (default: projected WGS84).
}
}

\value{
A \code{SpatialPointsDataFrame} containing the dataframe \code{occurrence}.
}

\seealso{
\code{\link{SpatialPointsDataFrame}}, \code{\link{checkOccurrence}}
}
\examples{
# load the occurrence dataframe
data(occurrence)

class(occurrence)
head(occurrence)

#convert to an SPDF
occ <- occurrence2SPDF(occurrence)

class(occ)
head(occ)
}