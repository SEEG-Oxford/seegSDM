\name{combinePreds}
\alias{combinePreds}

\title{
Combine multiple prediction rasters
}

\description{
Given a multilayer raster \code{preds} where each layer is a single bootstrap or ensemble prediction, calculate means and quantiles (specified by \code{quantiles}) for each pixel. If a snowfall cluster is running, \code{parallel = TRUE} sets the function to run in parallel. 
}
\usage{
combinePreds(preds, quantiles = c(0.025, 0.975), parallel = FALSE, maxn = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{preds}{
A \code{\link{RasterBrick}} or \code{\link{RasterStack}} object, where each layer is a single prediction.
}
  \item{quantiles}{
A vector of length two giving the quantiles of the prediciton distribution to return.
}
  \item{parallel}{
Whether to run the function in parallel. Requires a \code{\link{snowfall}} cluster to be running.
}
  \item{maxn}{
If \code{parallel = TRUE}, the maximum number of pixels in each batch of the parallel execution \code{maxn}. If \code{NULL}, \code{maxn} is set to fill all the cpus once.
}
}

\value{
A three layer \code{RasterBrick} or \code{RasterStack} object (whichever \code{preds} was) giving the mean and quantiles of the predictive distribution at each pixel.
}

\seealso{
\code{\link{RasterBrick}}, \code{\link{RasterStack}}, \code{\link{splitIdx}}, \code{\link{quantile}}, \code{\link{snowfall}}
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
preds_comb[[4]] <- preds_comb[[3]] - preds_comb[[2]]
names(preds_comb)[4] <- 'iq_range'

# plot them
plot(preds_comb)
}