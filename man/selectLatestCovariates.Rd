\name{selectLatestCovariates}
\alias{selectLatestCovariates}

\title{
  Select the lastest layer for each covariate in a nested list.
}
\description{
 From a named list containing synoptic and temporally resolved covariate layers, build a stack contain only the latest layer for each covariate (highest date).
}
\usage{
selectLatestCovariates(covariates)
}

\arguments{
  \item{covariates}{
    A named list of containing \code{RasterLayer} objects (or file paths), or sublists of \code{RasterLayer} (or file paths) for temporally resolved covariates. Synoptic covariates should have an entry with an approriate name and the \code{RasterLayer} (or file paths) as it's value. Annually resolved covariates should have an entry with an approriate name and a sub-list of files as it's value. The sub-list should name names in the format "YYYY", and values which are the \code{RasterLayer} (or file paths) of each covariate layer. Monthly resolved covariates should be structured in the same way as annual covariates, but use the format "YYYY-MM". Daily resolved covariates should be structured in the same way as annual covariates, but use the format "YYYY-MM-DD". 
}
}

\value{
  A \code{RasterStack} with a named layer for each covariate.
}

\seealso{
\code{\link{extractBatch}}
}
\examples{
# load occurrence data, admin units and covariates
data(covariates)
layers <- list(
  "cov1"=raster(covariates, layer="cov_a"),
  "cov2"=raster(covariates, layer="cov_b"),
  "cov3"=list(
    "2012"=raster(covariates, layer="cov_a"),
    "2013"=raster(covariates, layer="cov_b"),
    "2014"=raster(covariates, layer="cov_c")
  ),
  "cov4"=list(
    "2012"=raster(covariates, layer="cov_c"),
    "2013"=raster(covariates, layer="cov_b"),
    "2014"=raster(covariates, layer="cov_a")
  ),
  "cov5"=list(
    "2012-01"=raster(covariates, layer="cov_a"),
    "2012-02"=raster(covariates, layer="cov_b"),
    "2012-03"=raster(covariates, layer="cov_c")
  ),
  "cov6"=list(
    "2012-01"=raster(covariates, layer="cov_a"),
    "2012-02"=raster(covariates, layer="cov_b"),
    "2012-03"=raster(covariates, layer="cov_c")
  )
)

# extract all covariates,
stack <- selectLatestCovariates(layers)
}