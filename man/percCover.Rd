\name{percCover}
\alias{percCover}

\title{
Calculate the percentage cover of raster classes in a lower-resolution raster image.
}
\description{
Given discrete high resolution raster (\code{raster}) and a lower resolution template raster (\code{template}) and set of evaluation points, calculate the percentage cover of each class of \code{raster} in the cells of \code{template} identified by points.
}
\usage{
percCover(raster, template, points, codes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raster}{
A high-resolution discrete raster layer.
}
  \item{template}{
A lower resolution raster for which percentage cover of \code{raster} is to be calculated.
}
  \item{points}{
A set of points determining the pixels in \code{template} for which percentage cover is to be calculated. These are passed to the \code{x} argument of \code{\link{rasterize}}, so anything that function will accept should work here.
}
  \item{codes}{
A numeric vector giving the class labels for \code{raster}.
}
}
\value{
A dataframe giving the percentage cover in the required pixels.
}

\seealso{
\code{\link{rasterize}}
}
\examples{

# create a raster
r <- raster(nrows = 100, ncols = 100,
              xmn = 0, xmx = 10,
              ymn = 0, ymx = 10)
r <- setValues(r, round(runif(ncell(r)) * 10))
tmp <- aggregate(r, 10)
pts <- xyFromCell(tmp, sample(1:ncell(tmp), 3))
percCover(r, tmp, pts, 1:10)

}
