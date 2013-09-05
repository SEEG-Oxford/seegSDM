\name{checkOccurrences}
\alias{checkOccurrences}

\title{
  Quality control for input occurrence data
}

\description{
Given a dataframe \code{occurrence} of occurrence records and a raster giving the evidence \code{consensus}, run a series of quality control check on the data.
\enumerate{
  \item{
    The column names of \code{occurrence} are checked to ensure that the following required fields are present and of the correct class:
    \itemize{
      \item{'UniqueID' (integer)}{the record's unique identifying number}
      \item{'Admin' (integer, either 0, 1, 2, 3, or -9999)}{The admin level for polygons or -9999 for points}
      \item{'x' (numeric)}{the x coordinate (longitude) either of a point or a polygon centroid}
      \item{'y' (numeric)}{the y coordinate (latitude) either of a point or a polygon centroid}
      \item{'Area' (numeric)}{the area in square decimal degees of the polygon, or \code{NA} for a pixel}
    }
  }
  \item{Polygons with an 'Area' value greater than \code{area_threshold} square decimal degrees are removed.}
  \item{Records with coordinates which fall in cells with an evidence consensus value below \code{evidence_threshold} are removed.}
  \item{Coordinates are checked to make sure they don't fall outside the mask (i.e. in \code{consensus} cells which are \code{NA}). If any do fall in \code{NA}, they are moved to the nearest non-\code{NA} cell, provided it is less than \code{mas_distance} decimal degrees away. If this isn't possible the record is removed.}
}

If the occurrence data fails any checks outright, the function will stop and throw an error message. Otherwise the cleaned and corrected occurrence dataframe will be returned. \code{consensus} is assumed to conform to the mastergrid template and have a projected wgs84 coordinate reference. It may be worth checking this using \code{\link{checkRasters}} first.
}

\usage{
checkOccurrences(occurrences, consensus, consensus_threshold = -25,
                             area_threshold = 1, max_distance = 0.05, spatial = TRUE,
                             verbose = TRUE) 
}

\arguments{
  \item{occurrences}{a dataframe containing details of occurrence records}
  \item{consensus}{a \code{RasterLayer} object with projected WGS 1984 coordinate reference giving the regional evidence consensus scores (between -100 and 100). This should conform to the mastergrid template}
  \item{consensus_threshold}{a minimum evidence consensus value. Occurrence records will be removed if they fall in a cell of \code{consensus} with a value less than this}
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
# TO DO
}

