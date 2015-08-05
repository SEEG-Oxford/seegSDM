\name{extractBhatt}
\alias{extractBhatt}

\title{
Sample covariates according using a vector of three parameters
}
\description{
Generate pseudo-absence and pseudo-presence points according the procedure described by Bhatt et al. (2013) using the parameters: \code{na}, \code{np} and \code{mu} (specified by \code{pars}) and extract covariate values into a dataframe for modelling.

Where \code{n} is the number of occurrence points, \code{na * n} pseudo-absence points and \code{np * n} pseudo-presence points are generated at a maximum distance of \code{mu} decimal degrees from the occurrence points \code{occurrence}. The selection of pseudo-absence/presence points is spatially weighted based on the evidence consensus layer \code{consensus} (with values ranging from -100 to 100). Pseudo-absence points are more likely to be selected from cells with low consensus values (and never from cells with consensus 100). Pseudo-presence points are more likely to be from cells with high consensus values, and never from cells with values below \code{threshold}.

The covariate values for occurrence, pseudo-absence and pseudo-presence points are extracted from \code{covariates} and combined into a dataframe ready for modelling. Covariate values for polygon records are extracted using \code{\link{extractAdmin}} and summarized either by their mean or mode for continuous and discrete variables respectively. Discrete variables should by identified using the \code{factors} argument which will also coerce the covariates into factors in the resulting dataframe. If \code{points = TRUE} a list is returned containing the dataframe and the pseudo-presence and pseudo-absence points. If \code{return_points = FALSE} only the dataframe is returned.
}
\usage{
extractBhatt(pars, occurrence, covariates, consensus, admin, threshold = -25,
             factor = rep(FALSE, nlayers(covariates)), return_points = FALSE, ...)
}

\arguments{
  \item{pars}{
    Vector of length three giving parameters for selection of pseudo-presence and pseudo-absence points. The parameters are:
    \code{pars[1]}: the number of pseudo-absences to select for each occurrence point (na in Bhatt et al.)
    \code{pars[2]}: the number of pseudo-presences to select for each occurrence point (np in Bhatt et al.)
    \code{pars[3]}: the maximum distance (in decimal degrees) from occurrence points at which to select pseudo-presences and pseudo-absences (mu in Bhatt et al.)
  }
  \item{occurrence}{
    A \code{SpatialPointsDataFrame} object giving occurrence records and containing \emph{at least} the columns \code{Admin} (giving the admin level for polygons or -999 for points) and \code{GAUL} (giving the GAUL codes for polygons or \code{NA} for points).
  }
  \item{covariates}{
    A \code{Raster*} object giving covariates to use in the model.
  }
  \item{consensus}{
    A \code{RasterLayer} object giving the evidence consensus values (between -100 and 100).
  }
  \item{admin}{
  A \code{RasterBrick} or \code{RasterStack} object with four layers giving the GAUL codes for different admin
  levels. The layers must be in the order 0, 1, 2, 3, as they are in the example object \code{admin}.
  }
  \item{threshold}{
    A threshold evidence consensus value below which pseudo-presence points will not be selected.
  }
  \item{factor}{
    A logical vector stating whether each layer in \code{covariates} is discrete and should therefore be coerced into a factor.
  }
  \item{return_points}{
  Whether to return a list containing the dataframe and \code{SpatialPoints} objects of the pseudo-presence and pseudo-absence data (if \code{return_points = TRUE}) or just the dataframe (if \code{return_points = FALSE}).
  }
  \item{\dots}{
    Additional arguments to pass to \code{bgSample}. Note that \code{prob} is already set as \code{TRUE}, so will be ignored. 
  }
}

\value{
  If \code{return_points = TRUE}, a list with three elements:
    \item{data }{A dataframe with: column 1 giving 1s for occurrence and pseudo-presence records and 0s for pseudo-absence records; columns 2 and 3 giving the coordinates of these records; and the following columns containing the extracted covariates.}
    \item{pseudo_absence }{A \code{SpatialPoints} object of the pseudo-absence points selected. If \code{pars[1]} was 0 or lower this will be \code{NULL}}
    \item{pseudo_presence }{A \code{SpatialPoints} object of the pseudo-presence points selected. If \code{pars[2]} was 0 or lower this will be \code{NULL}}
    
  Otherwise just \code{data}.
}
\references{

Bhatt et al. (2013) The global distribution and burden of dengue. Nature \url{http://www.nature.com/nature/journal/v496/n7446/full/nature12060.html}
}

\seealso{
\code{\link{bgDistance}}, \code{\link{bgSample}}
}
\examples{
# load fake occurrence points, covariates and evidence consensus layer
data(occurrence)
data(covariates)
data(consensus)
data(admin)

# run checks, add GAUL codes and return a SpatialPoints object
occurrence <- checkOccurrence(occurrence, consensus, admin)

# run extractBhatt with:
# pseudo-absence ratio of 10:1
# pseudo-presence ratio of 0.1:1
# distance of 5 decimal degrees
# and specifying the last covariate as a discrete variable

lis <- extractBhatt(c(2, 1, 5),
                    occurrence,
                    covariates,
                    consensus,
                    admin,
                    return_points = TRUE,
                    factor = c(FALSE, FALSE, TRUE))

# look at what's produced
str(lis)

# specifically the dataframe
summary(lis$data)

# plot the true and pseudo points

# evidence consensus layer as a background
plot(consensus)
# add the pseudo-absences in light blue (note they aren't in high scoring
# consensus regions)
points(lis$pseudo_absence, pch = 16, col = 'light blue')
# and pseudo-presences in purple (not in the very low scoring regions) 
points(lis$pseudo_presence, pch = 16, col = 'purple')
# and add the occurrence points
points(occurrence, pch = 16)

}
