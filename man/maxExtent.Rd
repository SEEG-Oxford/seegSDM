\name{maxExtent}
\alias{maxExtent}

\title{
Create an extent object encompassing a list of corrdinate matrices.
}
\description{
Given a list (\code{list}) of matrices giving coordinates of points, calculate the smallest extent containing them all. Optionally, add a margin around the most extreme points using the \code{margin} argument.
}

\usage{
maxExtent(list, margin = c(0, 0, 0, 0))
}

\arguments{
  \item{list}{
  A list of matrices, each giving the latitudes (column 1) and longitudes (column 2) of a set of points.
}
  \item{margin}{
  An optional vector of length four giving a margin (in the map units) around the most extreme points.
}
}


\value{
An \code{extent} object encompassing the points in \code{list}.
}


\seealso{
\code{\link{extent}}
}
\examples{

# create coordinate matrices
one <- matrix(rnorm(3 * 2), ncol = 2)
two <- matrix(rnorm(5 * 2), ncol = 2)

e <- maxExtent(list(one, two), margin = c(1, 1, 1, 1))

}