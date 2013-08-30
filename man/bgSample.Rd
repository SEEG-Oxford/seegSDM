\name{bgSample}
\alias{bgSample}

\title{
Randomly sample pseudo-absence points from a raster
}
\description{
  Randomly sample \code{n} pseudo-absence points at pixels of \code{raster} (excluding \code{NA} pixels). If \code{prob = TRUE} the is assumed to be a bias grid and the random sampling is weighted toward pixels with the greatest value. If \code{spatial = TRUE} a \code{\link{SpatialPoints}} object is returned, else a matrix of coordinates.}

\usage{
bgSample(raster, n = 1000, prob = FALSE, replace = FALSE, spatial = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raster}{
The raster from which to sample at random.
}
  \item{n}{
The number of pseudo-absence points to sample.
}
  \item{prob}{
  Whether to weight samples by the values of the pixels in \code{raster}.
  }
  \item{replace}{
  Whether to sample with replacement.
  }
  \item{spatial}{
  Whether to return a \code{SpatialPoints} object (if \code{TRUE}, the default) or a matrix of coordinates (if \code{FALSE})
}
}
\value{
  Either a \code{SpatialPoints} object or a matrix of coordinates.
}

\seealso{
\code{\link{sample}}, \code{\link{bgDistance}}, \code{\link{biasGrid}}, \code{\link{SpatialPoints}}


}
\examples{
# load a test raster
r <- raster(system.file("external/test.grd", package="raster"))

# randomly sample 10000 background points without replacement
pts <- bgSample(r, n = 1000)

# plot the raster
plot(r)
# the pseudo-absences (with transparency)
plot(pts, add = TRUE, pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.2))


# this time with replacement 
pts <- bgSample(r, n = 1000, replace = TRUE)

plot(r)
plot(pts, add = TRUE, pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.2))


# now with replacement and biased by the values of the raster
pts <- bgSample(r, n = 1000, replace = TRUE, prob = TRUE)

plot(r)
plot(pts, add = TRUE, pch = 16, cex = 0.8, col = rgb(0, 0, 0, 0.2))
}
