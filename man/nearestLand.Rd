\name{nearestLand}
\alias{nearestLand}

\title{
Find nearest non-NA pixels to a set of points
}
\description{
Given a set of points (\code{points}), find the nearest cells in \code{raster} which are not \code{NA}, up to a maximum search radius of \code{max_distance}.
}
\usage{
nearestLand(points, raster, max_distance)
}

\arguments{
  \item{points}{
    A two-column matrix or dataframe giving x and y coordinates, or a \code{SpatialPoints*} object.
}
  \item{raster}{
    A \code{Raster*} object from which to find the nearest non-cells.
}
  \item{max_distance}{
    The maximum distance to search before giving up and returning \code{NA}. If \code{raster} is projected, this is in the map units, otherwise it is in metres.
}
}

\value{
  A two-column matrix giving the coordinates of the nearest cells, or \code{NA} values if there was no non-\code{NA} cell within \code{max_distance} of the points.
}

\seealso{
\code{\link{checkOccurrences}}
}

\examples{
# load the template raster
data(template)

# generate 100 random points on the raster
n <- 100
pts <- matrix(c(runif(n, -10, 10), runif(n, -5, 15)), ncol = 2)

# plot them
plot(template)
points(pts, pch = 16)

# find those which fall in NA areas on the raster
vals <- extract(template, pts)
outside_mask <- is.na(vals)
outside_pts <- pts[outside_mask, ]

# plot these in blue
points(outside_pts, pch = 16, col = 'blue')

# find the nearest land within 5 decimal degrees of these
land <- nearestLand(outside_pts, template, 2)

# count how many were moved
sum(!is.na(land[, 1]))

# and how many were too far out
sum(is.na(land[, 1]))

# plot the new points (for those which were reassigned) in green
points(land, pch = 16, col = 'dark green')

# plot where they moved to
arrows(outside_pts[, 1],
       outside_pts[, 2],
       land[, 1],
       land[, 2], length = 0.1)
}