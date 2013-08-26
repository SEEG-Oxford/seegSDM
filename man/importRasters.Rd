\name{importRasters}
\alias{importRasters}

\title{
Conveniently import raster files
}
\description{
Look in \code{path} and import all of the raster files with extension \code{ext} into a single object, determined by \code{as}. this function used \code{\link{mixedsort}} from the \code{gtools} package to pull in the files in alphanumeric order. 
}
\usage{
importRasters(path, as = brick, ext = ".grd")
}

\arguments{
  \item{path}{
  The file path in which to look for raster files.
}
  \item{as}{
  The function used to combine the rasters. \code{\link{brick}}, \code{\link{stack}} and \code{\link{list}} should all work.
}
  \item{ext}{
  The file extension for to identify rasters in \code{path}.
}
}

\value{
A single object, of the class determined by \code{as}.
}

\seealso{
\code{\link{mixedsort}}
}