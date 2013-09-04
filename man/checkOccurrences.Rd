\name{checkOccurrences}
\alias{checkOccurrences}

\title{
  Quality control for input occurrence data
}

\description{
Given a dataframe \code{occurrence} of occurrence records and a \code{RasterLayer} object \code{evidence} giving the evidence consensus map, this function:
\enumerate{
  \item{checks that the column names of \code{occurrence} match those expected (see details)}
  \item{checks that the data types of the columns in \code{occurrence} match those expected (see details)}
  \item{removes records of polygons with areas greater than \code{area_threshold}}
  \item{removes records ion areas with evidence consensus values below \code{evidence_threshold}}
  \item{checks for points outside the mask. If found, they will be moved to the nearest pixel inside the mask, provided it is no more than two pixels away. If this isn't possible the record is removed.}
}

If the occurrence data fails any checks outright (ie. stages 1:2), the function will stop and throw an error message. Otherwise the cleaned and corrected occurrence dataframe will be returned. \code{evidence} is assumed to conform to the mastergrid template. It may be worth checking this using \code{\link{checkRasters}} first.
  }

\usage{
checkOccurrences(occurrences, evidence, evidence_threshold = -100,
                  area_threshold = 1, verbose = TRUE)
}

\details{

TO DO:

expected fields and datatypes

evidence consensus for polygons based on point in field, not whole polygon

The input occurrence data

}

\arguments{
  \item{occurrences}{}
  \item{evidence}{}
  \item{evidence_threshold}{}
  \item{area_threshold}{}
  \item{verbose}{}
}

\value{
a list or something
}

\seealso{
\code{\link{checkRasters}}
}

\examples{
#
}

