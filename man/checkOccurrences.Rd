\name{checkOccurrences}
\alias{checkOccurrences}

\title{
  Quality control for input occurrence data
}

\description{
Given a dataframe \code{occurrence} of occurrence records and a \code{RasterLayer} object \code{evidence} giving the evidence consensus map, this function:
\begin{enumerate}
  \item checks that the column names of \code{occurrence} match those expected (see details)
  \item checks that the data types of the columns in \code{occurrence} match those expected (see details)
  \item removes records of polygons with areas greater than \code{area_threshold}
  \item removes records ion areas with evidence consensus values below \code{evidence_threshold}
  \item checks for points outside the mask. If found, they will be moved to the nearest pixel inside the mask, provided it is no more than two pixels away. If this isn't possible the record is removed.
\end{enumerate}
  If the occurrence data fails any checks outright (ie. stages 1:2), the function will stop and throw an error message.
  Otherwise the cleaned and corrected occurrence dataframe will be returned.
  
  \code{evidence} is assumed to conform to the mastergrid template. It may be worth checking this using \code{\link{checkRasters}} first.
  }

\usage{
checkOccurrences(occurrences, evidence, evidence_threshold = -100, area_threshold = 1, verbose = TRUE)
}

\details{

expected fields and datatypes


evidence consensus for polygons based on point in field, not whole polygon


The input occurrence data
In addition to the dataframe of covariates \code{x} a matrix of corresponding standard deviations, describing the measurement error, may be specified via the \code{error} argument.\cr\cr
If prior knowledge about the species' ecology is available (e.g. environmental thresholds beyond which the species cannot persist), this may be incorporated into the GRaF model by providing an R function via the \code{prior} argument.\cr\cr
The smoothness of the fitted random field is controlled by a vector of lengthscales which is used to parameterise a squared exponential covariance function, with one lengthscale for each covariate. \code{graf} centres and standardises covariates before fitting and the lengthscales therefore relate to these scaled covariates.
Lengthscales can be specified by the user via the the \code{l} argument, roughly approximated by \code{graf} by setting \code{l = NULL} or optimised by setting \code{opt.l = TRUE}.
}

\arguments{
  \item{occurrences}{}
  \item{evidence}{}
  \item{evidence_threshold}{}
  \item{area_threshold}{}
  \item{verbose}{}
}

\value{
A \code{graf} object which can be viewed and manipulated using \code{print}, \code{plot} and \code{predict} functions.
}

\author{
  Nick Golding, based on functions by Samir Bhatt
}

\references{
  Bhatt et al. (2013)  The global distribution and burden of dengue. Nature \url{http://dx.doi.org/10.1038/nature12060}
}

\seealso{
\code{\link{print.graf}}, \code{\link{plot.graf}}, \code{\link{predict.graf}}, \code{\link{plot3d}}, \code{\link{DIC}}, \code{\link{optim}}.
}

\examples{
# load Anguilla data from the dismo package
data(Anguilla_train)
}

