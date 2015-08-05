\name{rasterPointCount}
\alias{rasterPointCount}

\title{
Calculate the number of points falling within each pixel in a raster.
}
\description{
Given a rasterbrick and a two-column matrix of coordinates 'coords' (in the order: x, y), the function counts the number of points falling within each pixel in the rasterbrick and returns a two-column dataframe containing the pixel ID number and the count of points for that pixel.
}
\usage{
rasterPointCount(rasterbrick, coords)
}
\arguments{
  \item{rasterbrick}{
The raster from which to calculate the number of points within each pixel.
}
  \item{coords}{
A two-column matrix of coordinates in the order x then y.
}
}
\value{
A two-column dataframe with each row containing the pixel ID number and the count of points for that pixel. }
\seealso{
\code{\link{cellFromXY}}
}
\examples{
# load a test raster
r <- raster(system.file("external/test.grd", package="raster"))

# make some occurrence points
pts <- matrix(c(179000, 330000,
                181000, 333000,
                180000, 332000,
                180000, 332000),
                ncol = 2,
                byrow = TRUE)

# get count of points for each raster pixel
count_data <- rasterPointCount (r, pts)

}

