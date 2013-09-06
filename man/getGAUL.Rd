\name{getGAUL}
\alias{getGAUL}

\title{
Extract GAUL codes for polygon occurrence data
}
\description{
  Given a \code{SpatialPointsDataFrame} (\code{occurrence}) and a \code{RasterBrick} or \code{RasterStack} object containing raster layers of the four admin levels, extract GAUL codes for polygon records.
}
\usage{
getGAUL(occurrence, admin)
}

\arguments{
  \item{occurrence}{
  A \code{SpatialPointsDataFrame} of occurrence records containing \emph{at least} a column named \code{Admin} giving the admin level of the polygon. Admin should be an integer from 0 to 3, or -9999 for point data.
}
  \item{admin}{
  A \code{RasterBrick} or \code{RasterStack} object with four layers giving the GAUL codes for different admin levels. The layers must be in the order 0, 1, 2, 3, as they are in the example object \code{admin}.
}
}

\value{
A \code{SpatialPointsDataFrame} identical to \code{occurrence}, with an additional column \code{GAUL} giving the GAUL codes for the polygons and \code{NA} for points.
}

\seealso{
\code{\link{SpatialPointsDataFrame}}, \code{\link{checkOccurrence}}, \code{\link{occurrence2SPDF}}
}
\examples{
# load the occurrence and admin example objects
data(occurrence)
data(admin)

# convert occurrence into an SPDF
occurrence <- occurrence2SPDF(occurrence)

# look at occurrence
head(occurrence)

# note the order of admin
names(admin)

# extract GAUL codes
occurrence <- getGAUL(occurrence, admin)

# the GAUL colun has been added
head(occurrence)
tail(occurrence)
}