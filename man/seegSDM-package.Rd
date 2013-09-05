\name{seegSDM-package}
\alias{seegSDM-package}
\alias{seegSDM}
\docType{package}
\title{
Streamlined functions for species distribution modelling in the seeg research group.
}
\description{
The present focus is on ensemble boosted regression tree methods such as those used in Bhatt et al. (2013).
This will be developed to incorporate other statistical models, methods for selection of pseudo-absence data and other modelling choices.
}
\details{
\tabular{ll}{
Package: \tab seegSDM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-17\cr
License: \tab What license is it under?\cr
}
Functions are provided to assist with the seeg SDM workflow:
\enumerate{
  \item{quality control of occurrence and covariate data}
  \item{generation of pseudo-absence data}
  \item{fitting ensembles of models in parallel}
  \item{combining model ensembles to produce final predictions, uncertainty estimates and model summaries}
}
}
\author{
Nick Golding, with some functions based on code by Samir Bhatt
}
\references{
Bhatt et al. (2013)  The global distribution and burden of dengue. Nature \url{http://dx.doi.org/10.1038/nature12060}
}

\seealso{
Functions for manipulating rasters and shapefiles:

\code{\link{biasGrid}}, \code{\link{featureDensity}}, \code{\link{gaussWindow}}, \code{\link{osgb36}}, \code{\link{wgs84}}, \code{\link{en2os}}, \code{\link{extent2poly}}, \code{\link{importRasters}}, \code{\link{maxExtent}}, \code{\link{percCover}}, \code{\link{buildSP}}, \code{\link{nearestLand}}

Functions for generating pseudo-absences:

\code{\link{bgDistance}}, \code{\link{bgSample}}

Functions for running and assessing BRT ensembles:

\code{\link{runBRT}}, \code{\link{getEffectPlots}}, \code{\link{getRelInf}}, \code{\link{combinePreds}}

Functions for quality control:

\code{\link{checkOccurrences}}, \code{\link{checkRasters}}
 
Miscellaneous functions:

\code{\link{rmse}} \code{\link{sdWeighted}}, \code{\link{splitIdx}}, \code{\link{subsample}}

}