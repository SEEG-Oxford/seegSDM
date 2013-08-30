\name{bgDistance}
\alias{bgDistance}

\title{
Randomly sample pseudo-absence points from within a given distance of a set points
}
\description{
  Sample \code{n} background points from within a region within \code{distance} of \code{points} (a \code{\link{SpatialPoints}*} object)and in areas where \code{raster} is not \code{NA}.
  
The dots argument can be used to pass additional options to \code{\link{bgSample}}. To draw a biased sample (pseudo-absences more likely to be drawn from some areas within the distance buffer than others) set \code{raster} as the bias grid and add \code{prob = TRUE} via the dots argument (see examples for details).
}

\usage{
bgDistance(n, points, raster, distance, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  The number of pseudo-absence points to sample.
}
  \item{points}{
  A \code{SpatialPoints*} object from which to draw the buffer zone.
}
  \item{raster}{
  A template raster giving the regions from which to sample (non-NA areas). This can optionally be a bias grid (see examples).
}
  \item{distance}{
  The distance in map units defining the buffer zone around \code{points} from which to sample.
}
  \item{\dots}{
  Additional arguments to pass to \code{bgSample}.
}
}

\value{

A \code{SpatialPoints} object with \code{n} pseudo-absence points.

}

\seealso{
\code{\link{bgSample}}, \code{\link{buffer}}, \code{\link{biasGrid}}, \code{\link{SpatialPoints}}
}
\examples{
rm(list = ls())

# load a test raster
r <- raster(system.file("external/test.grd", package="raster"))

# find the extent
e <- extent(r)

# define some occurrence points
pts <- matrix(c(179000, 330000,
                181000, 333000,
                180000, 332000),
              ncol = 2,
              byrow = TRUE)

# build a SpatialPoints object with them
pts <- SpatialPoints(pts, r@crs)



# randomly sample 1000 background points from within 500m of presence points
# *without* replacement
pabs <- bgDistance(1000, pts, r, 500)

# plot the raster
plot(r)
# the pseudo-absences (with transparency)
plot(pabs, add = TRUE, pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.2))
# and the presences
plot(pts, add = TRUE, pch = 16, cex = 1.5, col = 'blue')



# this time sample with replacement
pabs <- bgDistance(1000, pts, r, 500, replace = TRUE)

plot(r)
plot(pabs, add = TRUE, pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.2))
plot(pts, add = TRUE, pch = 16, cex = 1.5, col = 'blue')



# now with replacement and biased by the values of the raster
pabs <- bgDistance(1000, pts, r ^ 2, 500, prob = TRUE, replace = TRUE)

plot(r)
plot(pabs, add = TRUE, pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.2))
plot(pts, add = TRUE, pch = 16, cex = 1.5, col = 'blue')

}
