\name{extractBatch}
\alias{extractBatch}

\title{
 Extract a batch of occurrence data. Takes account of synoptic or temporally resolved covariates, as well as, point and admin data.
}
\description{
Given a \code{SpatialPointsDataFrame} object \code{batch} with \code{Admin}, \code{GAUL} and \code{Date} columns and a named list of covariate file paths extract the values of \code{covariates}.
}
\usage{
extractBatch(batch, covariates, factor, admin, admin_mode="average", load_stack=stack)
}

\arguments{
  \item{batch}{
A \code{SpatialPointsDataFrame} of occurrence records containing \emph{at least} a column named \code{Admin} giving the admin level of the occurrence (an integer from 0 to 3, or -999 for point data), a column named \code{GAUL} giving the GAUL codes for polygons occurrences and a column named \code{Date} giving the occurrence date of the occurrence as a Date object.
}
  \item{covariates}{
    A named list of containing \code{RasterLayer} objects (or file paths), or sublists of \code{RasterLayer} (or file paths) for temporally resolved covariates. Synoptic covariates should have an entry with an approriate name and the \code{RasterLayer} (or file paths) as it's value. Annually resolved covariates should have an entry with an approriate name and a sub-list of files as it's value. The sub-list should name names in the format "YYYY", and values which are the \code{RasterLayer} (or file paths) of each covariate layer. Monthly resolved covariates should be structured in the same way as annual covariates, but use the format "YYYY-MM". Daily resolved covariates should be structured in the same way as annual covariates, but use the format "YYYY-MM-DD". 
}
  \item{factor}{
    A named list stating whether each entry in \code{covariates} is discrete and should therefore be coerced into a factor.
  }
  
  \item{admin}{
  A \code{RasterBrick} or \code{RasterStack} object with four layers giving the GAUL codes for different admin levels. The layers must be in the order 0, 1, 2, 3, as they are in the example object \code{admin}.
}
  \item{admin_mode}{
  A string indicating the type of extraction to perform for polygon occurrences. Acceptable values are: "average" (returns the mean value for non-discrete covariates, and the mode for discrete covariates), "random" (returns the value of a random pixel) and "latlong" (returns the value at the lat/long in the \code{SpatialPointsDataFrame}, normally the centroid).
}
  \item{load_stack}{
  (Optional) A function which creates a stack from list of file paths. This can be used to override the crs and/or extent read from the files. If "NA" is specified, raster::stack will be used. 
}
}

\value{
  A matrix giving the covariate values for each of the points/admins in \code{batch}.
}

\seealso{
\code{\link{extractAdmin}}
}
\examples{
# load occurrence data, admin units and covariates
# normally these would be individual geotiffs on disk
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
factor <- list(
  "cov1"=FALSE,
  "cov2"=TRUE,
  "cov3"=FALSE, 
  "cov4"=TRUE,
  "cov5"=FALSE,
  "cov6"=TRUE
)

data(admin)

data(occurrence)
occurrence <- occurrence2SPDF(occurrence)
occurrence <- getGAUL(occurrence, admin)
occurrence$Date <- as.Date(sample(c(
  "2011-10-01", 
  "2012-01-02", 
  "2012-02-02", 
  "2012-03-06", 
  "2012-04-06",
  "2015-03-06"), length(occurrence), replace=TRUE))
occurrence$PA <- sample(c(0, 1), length(occurrence), replace=TRUE)

# extract all covariates,
vals <- extractBatch(occurrence, layers, factor, admin, admin_mode="average")
head(vals)
}