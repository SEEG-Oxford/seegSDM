\name{Eastings/Northings}
\alias{en2os}
\alias{os2en}
\alias{en2poly}


\title{
Converting Ordnance Survey grid references
}

\description{
Convert between alphanumeric Ordnance Survey grid references and all-numeric eastings/northings for the Ordnance Survey of Great Britain coordinate system (\code{os2en} & \code{en2os})and convert Eastings/Northings into polygons (\code{en2poly}).
}

\usage{
en2os(en)
os2en(os)
en2poly(coords)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{en}{
A two-column matrix of Ordnance Survey eastings (column 1) and northings (column 2) coordinates of arbitrary precision.
}
  \item{os}{
A character vector of alphanumeric Ordnance Survey grid references of arbitrary precision, optionally with DINTY tetrad letters (see references).
}
  \item{coords}{
A matrix with four columns giving (in this order) eastings of the lower-left corners of the grid squares, northings of the lower-left corners of the grid squares, eastings of the upper-right corners of the grid square, northings of the upper-right corners of the grid squares.
}
}

\value{
\code{en2os}: a character vector of alphanumeric OS grid references.

\code{os2en}: a dataframe giving the lower-left and upper right corners and centre of the OS grid square in eastings/northings and the resolution (square side length) in metres.

\code{en2poly}: a \code{\link{SpatialPolygons}} object with a polygon for each grid square.
}

\references{
Ordnance Survey \url{http://en.wikipedia.org/wiki/Ordnance_Survey}
'DINTY' tetrads \url{http://en.wikipedia.org/wiki/Tetrad_(area)}
}

\seealso{
\code{\link{osgb36}}
}

\examples{
# generate some OS grid references of differing precisions (note the third has a DINTY tetrad reference)
os <- c('SP320254', 'SP3209225422', 'SP320254H')
# convert them to eastings & Northings
en <- os2en(os)
en
# and back to grid references
os2 <- en2os(en[, 1:2])
# note that the eastings-northings sytem cannot handle dinty tetrads, so the information is lost
os2
# plot the polygons for the original OS coordinates
polys <- en2poly(en[, 1:4])
plot(polys, axes = TRUE)
}
