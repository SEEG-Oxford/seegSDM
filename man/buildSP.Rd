\name{buildSP}
\alias{buildSP}

\title{
Combine multiple shapefiles into a \code{SpatialPolygons} object
}
\description{
Given a set of filenames for shapefiles, load them all into a single \code{\link{SpatialPolygons}} object.
}
\usage{
buildSP(files)
}

\arguments{
  \item{files}{
A vector of filenames corresponding to shapefiles (ending in \code{.shp})
}
}

\value{
A \code{SpatialPolygons} object.
}

\seealso{
\code{\link{SpatialPolygons}}, \code{\link{shapefile}}
}
