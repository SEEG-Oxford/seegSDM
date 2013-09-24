\name{tempStand}
\alias{tempStand}

\title{
Temporal standardisation of occurrence data
}

\description{
Expand a dataframe of occurrence by year range so that there is a separate record for each year. Also check for duplicated year / location combinations.

\code{occurrence} should be a dataframe with columns named \code{Start_Year}, \code{End_Year}, \code{Admin}, \code{Latitude} and \code{Longitude} giving the start and end years of the date range, the admin level (0 to 3, or -999 for points) and the latitude and longitude of points or centroids. A \code{SourceID} column should also be present, giving a unique identifier for each of the input records. If a column of this name is not present one will be added. A \code{UniqueID} column is also added, giving a unique identifier for each record after the temporal standardisation.

GAUL codes for each record are required to check for duplicates. If these are provided as a column named \code{GAUL} in \code{occurrence}, then this is used, otherwise \code{\link{getGAUL}} is used to extract the GAUL codes from \code{admin} (a \code{RasterBrick} or \code{RasterStack} object giving the GAUL codes at different admin levels) and add them to \code{occurrence}.

\code{tempStand} returns a list giving the expanded dataframe and row indices for the duplicated points and polygons.
}
\usage{
tempStand(occurrence, admin, verbose = TRUE)
}

\arguments{
  \item{occurrence}{
  A dataframe of occurrence data containing the columns \code{Start_Year}, \code{End_Year}, \code{Admin}, \code{Latitude} and \code{Longitude} and preferably \code{SourceID} and \code{GAUL}.
}
  \item{admin}{
A \code{RasterBrick} or \code{RasterStack} object with four layers giving the GAUL codes for different admin levels. The layers must be in the order 0, 1, 2, 3, as they are in the example object \code{admin}.
}
  \item{verbose}{
  Whether to print information to the console.
}
}

\value{
A list with three elements:
  \item{occurrence }{An expanded dataframe of occurrence data with individual records for each year. The columns \code{Start_Year} and \code{End_Year} are removed and columns \code{Year} and \code{UniqueID} are added (as well as \code{SourceID} and \code{GAUL} if they weren't already present). Note that this dataframe still contains the duplicate records, which the user should deal with before continuing.
}

\item{duplicated_polygons }{A vector giving the row numbers of the output element \code{occurrence} which contain duplicated polygon / year combinations. Note that only the first identical element is not listed.}

\item{duplicated_points }{A vector giving the row numbers of the output element \code{occurrence} which contain duplicated point / year combinations. Points are considered to be duplicates if they fall in the same pixel of \code{admin}. Note that the first of the identical elements is not listed.}
}
\seealso{
\code{\link{getGAUL}}, \code{\link{duplicated}}
}
\examples{
# TO DO
}