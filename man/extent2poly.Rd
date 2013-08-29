\name{extent2poly}
\alias{extent2poly}

\title{
Convert an extent object into a polygon
}
\description{
Convert an \code{\link{extent}} object into a \code{\link{SpatialPolygons}} object. The coordinate reference is not recorded in the extent object so must aso be provided
}
\usage{
extent2poly(extent, proj4string)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{extent}{
    An \code{extent} object.
}
  \item{proj4string}{
    A valid coordinate reference (passed to the \code{proj4string} argument in \code{SpatialPolygons})
}
}

\value{
A rectangular \code{SpatialPolygons} object matching \code{extent}.
}

\seealso{
\code{\link{extent}}, \code{\link{SpatialPolygons}}
}

\examples{
ext <- extent(c(0, 10, 0, 5))
poly <- extent2poly(ext, wgs84())
plot(poly)
}