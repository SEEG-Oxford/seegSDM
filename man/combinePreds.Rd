\name{combinePreds}
\alias{combinePreds}

\title{
Combine multiple prediction rasters
}

\description{
Given a multilayer raster \code{preds} where each layer is a single bootstrap or ensemble prediction, calculate mean, median and and quantile (specified by \code{quantiles}) rasters.
}
\usage{
combinePreds(preds, quantiles = c(0.025, 0.975),
              parallel = FALSE, ncore = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{preds}{
A \code{RasterBrick} or \code{RasterStack} object, where each layer is a single prediction.
}
  \item{quantiles}{
A vector of length two giving the quantiles of the prediciton distribution to return.
}
  \item{parallel}{
Whether to run the function in parallel. The user can specify the number of cores to run on via the \code{ncore} argument, or if this is \code{NULL} and a \code{\link{snowfall}} cluster is running the existing cluster will be used. If \code{parallel = TRUE} but \code{ncore} is not provided and a \code{snowfall} cluster isn't running, only one core will be used.
}
  \item{ncore}{
If \code{parallel = TRUE} the number of cores to use, see \code{parallel} for details.
}
}

\value{
A \code{RasterBrick} or \code{RasterStack} object (whichever \code{preds} was) giving the mean, median and as many quantiles as were requested (default 2) of the predictive distribution at each pixel.
}

\seealso{
\code{RasterBrick}, \code{RasterStack}, \code{\link{quantile}}
}
\examples{
# load a test raster
r <- raster(system.file("external/test.grd", package="raster"))
# generate multiple layers with random perturbations
preds <- lapply(1:10, function(i, raster) {
  raster[] <- raster[] + rnorm(ncell(raster), 0, 10 * sqrt(raster[] + 1))
  raster
}, r)

# turn them into a brick
preds <- brick(preds)

# combine them, returning the mean and the 10% and 90% quantiles
preds_comb <- combinePreds(preds, quantiles = c(0.1, 0.9))

# calculate the inter-quantile range
preds_comb[[5]] <- preds_comb[[4]] - preds_comb[[3]]
names(preds_comb)[5] <- 'iq_range'

# plot them
plot(preds_comb)
}