\name{masterMask}
\alias{masterMask}

\title{
Create master mask from raster stack
}
\description{
Given a stack of rasters. Loop through the rasters masking one layer by every other layer to create a master mask.
}
\usage{
masterMask(rasters)
}

\arguments{
  \item{rasters}{
A raster stack.
}
}

\value{
Raster object.
}

\seealso{
\code{\link{mask}}, ~~~
}
\examples{
# load a raster stack with 3 layers
rs <- stack(system.file("external/rlogo.grd", package="raster"))

# create a master mask
master_mask <- masterMask(rs)

function (rasters) 
{
    master <- rasters[[1]]
    for (i in 1:nlayers(rasters)) {
        master <- mask(master, rasters[[i]])
    }
    return(master)
  }
}


