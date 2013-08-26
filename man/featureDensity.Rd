\name{featureDensity}
\alias{featureDensity}

\title{
Create a raster of the density of points and polygons
}
\description{
Given a feature object (\code{feature}) and a template raster (\code{raster}), create a raster giving the density of features at each cell. The density of polygons can be split accross multiple cells be providing a vector of weights (perhaps calculated using \code{\link{areaDensity}}) via the \code{weights} argument. This is essentially a wrapper to \code{\link{rasterize}} in the \code{raster} package.
}
\usage{
featureDensity(feature, raster, weights = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{feature}{
  A feature object, anything that is accepted as the \code{x} argument in \code{rasterize}.
}
  \item{raster}{
  A template raster for the resulting density raster.
}
  \item{weights}{
  An optional vector of weights, for the elements of \code{feature} or a string giving a column of \code{feature} containing weights. This is passed to the \code{field} argument in \code{rasterize}, so read that helpfile for more details.
}
  \item{\dots}{
Other arguments to pass to \code{rasterize}.
}
}

\value{
A raster object with the same dimension and resolution as \code{raster} giving the density of features.
}

\seealso{
\code{\link{rasterize}}, \code{\link{areaDensity}}
}
\examples{
#
# create a raster
r <- raster(nrows = 10, ncols = 10,
              xmn = 0, xmx = 10,
              ymn = 0, ymx = 10)
r <- setValues(r, round(runif(ncell(r)) * 10))

# and some points
pts <- xyFromCell(r,
                  sample(1:ncell(r), 500, replace = TRUE))

# some fake weights
wts <- runif(nrow(pts))
                  
# calculate their density
density <- featureDensity(pts, r, weights = wts)
}
