\name{biasGrid}
\alias{biasGrid}

\title{
  Create a bias grid from polygons using a gaussian moving window smoother
}

\description{
Given a \code{\link{SpatialPolygons}*} object representing observed occurrence of a reference disease group (\code{polygons}) and a template raster (\code{raster}), use a Gaussian moving window smoother with standard deviation \code{sigma} to generate a smoothed raster of the expected density of observations. See \code{\link{gaussWindow}} for details of the Gaussian smoother. The resulting bias grid can be used to select pseudo-absence points in order to account for observation bias in the occurrence dataset.
}

\usage{
biasGrid(polygons, raster, sigma = 30)
}

\arguments{
  \item{polygons}{
  A \code{SpatialPolygons*} object.
}
  \item{raster}{
  A template raster for the resulting bias grid.
}
  \item{sigma}{
  The standard deviation, in map units, of the Gaussian smoother.
}
}

\value{
A raster of the expected density of occurrence records.
}

\seealso{
\code{\link{gaussWindow}}, \code{\link{focal}}
}
\examples{
# load a test raster
raster <- raster(system.file("external/test.grd", package="raster"))

# define some occurrence points
pts <- matrix(c(179000, 330000,
                181000, 333000,
                180000, 332000),
              ncol = 2,
              byrow = TRUE)

# build a SpatialPoints object with them
pts <- SpatialPoints(pts, raster@crs)

# buffer them to create polygons
polygons <- gBuffer(pts, width = 50, byid = TRUE)

# run biasGrid with a fairly large sigma (500 metres)
bias <- biasGrid(polygons, raster, sigma = 500)

# and plot the results
plot(bias)
plot(polygons, add = TRUE)
}


