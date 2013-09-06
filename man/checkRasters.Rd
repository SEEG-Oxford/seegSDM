\name{checkRasters}
\alias{checkRasters}

\title{
Check that rasters confrom to a specified template
}
\description{
Check that the raster(s) \code{rasters} conforms to a raster template, specified by \code{template}. By default the extent, projection and number of cells are checked. If \code{celllbycell = TRUE} cell values in \code{rasters} and \code{template} are individually compared to ensure that all \code{NA} values in \code{template} are \code{NA} in \code{rasters} and vice versa - but this can be very slow with big rasters.

If any of the tests are failed \code{checkRasters} throws an error message (which should hopefully be informative) and stops, if they are passed \code{rasters} is returned. The function should therefore be used with care in scripts since calls of the form \code{ras <- checkRasters(ras, template)} could throw an error but allow the script to continue. Consider wrapping the workflow up in a function, using \code{\link{tryCatch}} or assigning to a new object, as detailed in the examples.
}
\usage{
checkRasters(rasters, template, cellbycell = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rasters}{
A \code{RasterLayer}, \code{RasterBrick} or\code{RasterStack} object containing the raster  or rasters to check.
}
  \item{template}{
A \code{RasterLayer} object giving the template raster to compare \code{rasters} to.
}
  \item{cellbycell}{
Whether to individually check cell values between \code{rasters} and \code{template} to ensure that \code{NA} values (masked values) line up. This can be very slow for large rasters!
}
}
\value{
If all checks are passed \code{rasters} is returned, otherwise \code{checkRasters} throws an error.
}

\seealso{
\code{\link{checkOccurrence}}, \code{\link{tryCatch}}
}
\examples{
## Not run:
# (not running example when checking package as
#   error messages would cause it to fail)

# load the example template raster
data(template)

# make a load of copys
raster_proj <- raster_extent <- raster_ncell <- raster_values <- template

# and break them!
# wrong projection
projection(raster_proj) <- wgs84()

# wrong extent
extent(raster_extent) <- extent(as.vector(extent(raster_extent)) / 2)

# wrong number of cells
raster_ncell <- aggregate(raster_ncell, 2)

# and missing cell values
raster_values[sample(1:nrow(template), 20), sample(1:ncol(template), 20)] <- NA
\dontrun{
# now test them
ans <- checkRasters(raster_proj, template)
ans <- checkRasters(raster_extent, template)
ans <- checkRasters(raster_ncell, template)
ans <- checkRasters(raster_values, template, cellbycell = TRUE)

# so nothing has been assigned to ans
ans

# note that using checkRasters to assign back to the same object is a bad idea
# in a script:
raster_proj <- checkRasters(raster_proj, template)
# since it returns an error message, but the script continues with the original
# dodgy raster as raster_proj wasn't overwritten by anything

# instead you could wrap up the whole analysis in a function like this:
run <- function (raster, template) {
  raster <- checkRasters(raster, template)
  plot(raster)
}

# since the function will fail with the error message and stop
run(raster_proj, template)

# (this is what it does with no error)
run(template, template)

# alternatively you could use tryCatch to overwrite with NULL
# if there's an error
raster_proj <- tryCatch(checkRasters(raster_proj, template),
                        error = function(error) NULL)
raster_proj

# An easier solution would be to assign the answer to a new object,
# causing the rest of the script to fail if the checks aren't completed
raster_clean <- checkRasters(raster_extent, template)

plot(raster_clean)
# doesn't exist so throws an error
}
## End(Not run)
}