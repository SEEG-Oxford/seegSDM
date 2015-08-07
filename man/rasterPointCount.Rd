\name{rasterPointCount}
\alias{rasterPointCount}

\title{
Calculate the number of points falling within each pixel in a raster.
}
\description{
Given a rasterbrick and a two-column matrix of coordinates 'coords' (in the order: x, y), the function counts the number of points falling within each pixel in the rasterbrick and returns a three-column dataframe containing the x and y coordinates for each pixel and the count of occurrence points for that pixel. Pixels containing no points are assigned the value 0. 

Option to include pseudo-absence points. A second two-column matrix of coordinates 'absence' (in the order: x, y)   

If \code{extract = TRUE} raster values for each pixel are extracted and also returned in the dataframe.

If any coordinates for occurrence and pseudo-absence points fall within the same pixel, the function returns an error.

}
\usage{
rasterPointCount(rasterbrick, coords, absence=NULL, extract=FALSE)
}
\arguments{
  \item{rasterbrick}{
The raster from which to calculate the number of points within each pixel.
}
  \item{coords}{
A two-column matrix of coordinates in the order x then y. Occurrence points.
}
 \item{absence}{
A two-column matrix of coordinates in the order x then y. Pseudo-absence points.
}
}
\value{
A three-column dataframe containing x and y coordinates for each pixel and the count of occurrence points for that pixel. Pixels containing no points are assigned the value 0. 

If pseudo-absence points are included, the function returns a three-column dataframe containing x and y coordinates for each pixel and the count of occurrence points for that pixel. Pixels containing pseudo-absence points are assigned the value 0. Pixels containing no points are assigned the value NA. 
}
\seealso{
\code{\link{cellFromXY}}
}
\examples{
# load a test raster
r <- raster(system.file("external/test.grd", package="raster"))

# make some occurrence points
occ <- matrix(c(179000, 330000,
                181000, 333000,
                180000, 332000,
                180000, 332000),
                ncol = 2,
                byrow = TRUE)

# make some pseudo-absence points
absence <- matrix(c(179000, 331000,
                    179000, 330500,
                    180000, 331000),
                    ncol = 2,
                    byrow = TRUE)

# get count of points for each raster pixel
count_data <- rasterPointCount (r, occ, absence)

}

