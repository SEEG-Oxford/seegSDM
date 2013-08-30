\name{gaussWindow}
\alias{gaussWindow}

\title{
Generate a Gaussian window for smoothing rasters.
}
\description{
Create a window (essentially just a matrix) of dimension\code{n} with Gaussian weights (using dispersion parameter \code{sigma}) for use in moving window smoothing (e.g. \code{\link{focal}} from the \code{raster} package or \code{\link{biasGrid}} from this package).
}
\usage{
gaussWindow(n, sigma = n/5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
The dimension of the window, i.e. the window has dimension \code{n * n}. \code{n} must be an odd number (in order to have the focal cell in the center) so if an even number is specified \code{n + 1} is used instead, without a warning.
}
  \item{sigma}{
The dispersion of the Gaussian.
}
}
\details{
The density of the Gaussian is calculated as:

  \code{exp(-(r ^ 2) / (2 * sigma ^ 2))}
  
where r is the euclidean distance (in number of cells) from the focal (center) cell.
}
\value{
A square matrix of dimension \code{n} (or \code{n + 1}, see above) giving the weights for moving window smoothing.
}

\seealso{
\code{\link{focal}}, \code{\link{biasGrid}}
}
\examples{
# create a raster
r <- raster(nrows = 10, ncols = 10,
              xmn = 0, xmx = 10,
              ymn = 0, ymx = 10)
r <- setValues(r, round(runif(ncell(r)) * 10))

# set up the Gaussian window, using the default value of sigma
window <- gaussWindow(n = 5)
window

# run the moving window smoother, padding and ignoring NAs
r2 <- focal(r, window, pad = TRUE, na.rm= TRUE)

# plot the before and after
par(mfrow = c(1, 2))
plot(r)
plot(r2)

}