\name{wgs84}
\alias{wgs84}

\title{
World Geodetic System 1984 coordinate reference
}
\description{
A helper function returning a \code{CRS} object giving the WGS84 coordinate reference, either projected or unprojected.
}
\usage{
wgs84(projected = FALSE)
}

\arguments{
  \item{projected}{
  Whether to return the projected WGS84 reference, defaults to the unprojected.
}
}

\value{
An object of class \code{CRS}
}

\seealso{
\code{\link{CRS-class}}
}

\examples{
wgs84()
# CRS arguments: +init=epsg:4326
wgs84(TRUE)
# CRS arguments: +init=epsg:3395

# create a raster
r <- raster(nrows = 10, ncols = 10,
              xmn = 0, xmx = 10)
# by default wgs84
projection(r)
# [1] "+proj=longlat +datum=WGS84"

# change it to projected (not that you'd want to do this i.r.l.)
projection(r) <- wgs84(TRUE)
projection(r)
#[1] "+init=epsg:3395"

}
