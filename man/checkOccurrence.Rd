\name{checkOccurrence}
\alias{checkOccurrence}

\title{
  Quality control for input occurrence data
}

\description{
Given a dataframe \code{occurrence} of occurrence records and rasters giving the evidence consensus (\code{consensus}) and admin units (\code{admin}), run a series of quality control check on the data. The checks are:
\enumerate{
  \item{
    The column names of \code{occurrence} are checked to ensure that the following required fields are present and of the correct class:
    \itemize{
      \item{}{'UniqueID' (integer): the record's unique identifying number}
      \item{}{'Admin' (integer, either 0, 1, 2, 3, or -999): the admin level for polygons or -999 for points}
      \item{}{'Year' (integer): the year of the occurrence record}
      \item{}{'x' (numeric): the x coordinate (longitude) either of a point or a polygon centroid}
      \item{}{'y' (numeric): the y coordinate (latitude) either of a point or a polygon centroid}
      \item{}{'Area' (numeric): the area in square decimal degees of the polygon, or \code{NA} for a pixel}
    }
    Any error is thrown if any of these are not presnt or of the wrong class.
  }
  \item{Polygon records are checked for duplicated polygon/year combinations and an error is thrown if any are present.}
  \item{Polygons with an 'Area' value greater than \code{area_threshold} square decimal degrees are removed.}
  \item{Records with coordinates which fall in cells with an evidence consensus value below \code{evidence_threshold} are removed.}
  \item{Coordinates are checked to make sure they don't fall outside the mask (i.e. in \code{consensus} cells which are \code{NA}). If any do fall in \code{NA}, they are moved to the nearest non-\code{NA} cell, provided it is less than \code{max_distance} decimal degrees away. If this isn't possible the record is removed.}
  
  The GAUL code for any polygons
}

If the occurrence data fails any checks outright, the function will stop and throw an error message. Otherwise the cleaned and corrected occurrence dataframe will be returned. \code{consensus} is assumed to conform to the mastergrid template and have a projected wgs84 coordinate reference. It may be worth checking this using \code{\link{checkRasters}} first.
}

\usage{
checkOccurrence(occurrence, consensus, admin, consensus_threshold = -25,
                             area_threshold = 1, max_distance = 0.05,
                             spatial = TRUE, verbose = TRUE) 
}

\arguments{
  \item{occurrence}{A dataframe containing details of occurrence records}
  \item{consensus}{A \code{RasterLayer} object with projected WGS 1984 coordinate reference giving the regional evidence consensus scores (between -100 and 100). This should conform to the mastergrid template}
  \item{admin}{A \code{RasterBrick} or \code{RasterStack} object with four layers giving the GAUL codes for different admin levels. The layers must be in the order 0, 1, 2, 3, as they are in the example object \code{admin}.}
  \item{consensus_threshold}{A minimum evidence consensus value. Occurrence records will be removed if they fall in a cell of \code{consensus} with a value less than this}
  \item{area_threshold}{The minimum area (in square decimal degrees) of polygons for inclusion. Polygon records with a lower value will be removed.}
  \item{max_distance}{The maximum distance (in decimal degrees) to search for a non-\code{NA} cell when attempting to reassign coordinates which fall in \code{NA} cells of \code{consensus}}
  \item{spatial}{Whether to return a \code{SpatialPointsDataFrame} (if \code{TRUE}) or just a dataframe (if \code{FALSE})}
  \item{verbose}{Whether to print information on non-critical issues to the console.}
}

\value{
Provided none of the checks were failed outright, either a dataframe (if \code{spatial = FALSE}) or a \code{SpatialPointsDataFrame} (if \code{spatial = TRUE}) containing the cleaned and checked occurrence data.
}

\seealso{
\code{\link{checkRasters}}, \code{\link{nearestLand}}
}

\examples{
# load occurrence data, consensus and admin objects
data(occurrence)
data(consensus)
data(admin)

head(occurrence)
class(occurrence)

# run checkOccurrence
occ <- checkOccurrence(occurrence, consensus, admin)

head(occ)
class(occ)
}

