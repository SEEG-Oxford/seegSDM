\name{extractAdmin}
\alias{extractAdmin}

\title{
 Extract covariates by admin units
}
\description{
Given a \code{SpatialPointsDataFrame} object \code{occurrence} with \code{Admin} and \code{GAUL} columns and a \code{RasterBrick} or \code{RasterStack} object \code{admin} giving the admin layers,  extract and summarize the values of \code{covariates} using the function \code{fun} for each polygon.
}
\usage{
extractAdmin(occurrence, covariates, admin, fun = mean)
}

\arguments{
  \item{occurrence}{
A \code{SpatialPointsDataFrame} of occurrence records containing \emph{at least} a column named \code{Admin} giving the admin level of the polygon (an integer from 0 to 3, or -999 for point data) and a column named \code{GAUL} giving the GAUL codes for these polygons. This can be added using \code{\link{getGAUL}}.
}
  \item{covariates}{
    A \code{Raster*} object giving covariates to use in the model.
}
  \item{admin}{
  A \code{RasterBrick} or \code{RasterStack} object with four layers giving the GAUL codes for different admin levels. The layers must be in the order 0, 1, 2, 3, as they are in the example object \code{admin}.
}
  \item{fun}{
  A function to summarize the covariate values accross each polygon, \code{mean} might be a good option for continuous covariates, or \code{modal} for discrete covariates.
}
}

\value{
  A matrix giving the summarized covariate values for each of the \emph{polygons} in \code{occurrence} (i.e. point records will be missing).
}

\seealso{
\code{\link{getGAUL}}
}
\examples{
# load occurrence data, admin units and covariates
data(occurrence)
data(admin)
data(covariates)

# convert occurrence to an SPDF
occurrence <- occurrence2SPDF(occurrence)
# add on the GAUL codes
occurrence <- getGAUL(occurrence, admin)

# extract all covariates, taking the mean across each polygon
vals <- extractAdmin(occurrence, covariates, admin, fun = mean)
head(vals)

# note that point records are omitted
nrow(occurrence)
nrow(vals)

# it would be better to use different summary functions since cov_c is a factor
cov_ab <- extractAdmin(occurrence, covariates[[1:2]], admin, fun = mean)
cov_c <- extractAdmin(occurrence, covariates[[3]], admin, fun = modal)

head(cbind(cov_ab, cov_c))
}