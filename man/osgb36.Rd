\name{osgb36}
\alias{osgb36}

\title{
British Ordnance Survey 1936 coordinate reference
}
\description{
A helper function returning a \code{CRS} object giving the OSGB36 coordinate reference.
}
\usage{
osgb36()
}

\value{
An object of class \code{CRS}
}

\seealso{
\code{\link{CRS-class}}
}

\examples{
osgb36()
# CRS arguments:
#  +proj=tmerc +lat_0=49 +lon_0=-2 +k=0.999601272 +x_0=400000 +y_0=-100000
# +datum=OSGB36 +units=m +no_defs +ellps=airy
# +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894 

# create a raster
r <- raster(nrows = 10, ncols = 10,
              xmn = 0, xmx = 10)
# by default wgs84
projection(r)
# [1] "+proj=longlat +datum=WGS84"

# change it to OSGB36
projection(r) <- osgb36()
projection(r)
# [1] "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.999601272 +x_0=400000\n  +y_0=-100000
# +datum=OSGB36 +units=m +no_defs +ellps=airy
# +towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894"
}
