\name{seegSDM synthetic data}
\alias{occurrence}
\alias{evaluation}
\alias{covariates}
\alias{consensus}
\alias{template}
\alias{admin}

\docType{data}
\title{
Synthetic data for demonstrating the seegSDM package
}

\description{
A set of artificial data for demonstrating and testing the \code{seegSDM} R package.
}

\usage{
data(occurrence)
data(evaluation)
data(covariates)
data(consensus)
data(template)
data(admin)
}

\format{
  \code{occurrence}: a dataframe giving occurrence data for a hypothetical disease (columns \code{UniqueID}, \code{Admin}, \code{Year}, \code{x}, \code{y} and \code{Area} - as required by \code{\link{checkOccurrence}}).
  
  \code{evaluation}: a matrix with three columns, giving (projected) lat longs of sites which have been surveyed for the hypothetical disease (\code{x} and \code{y}) and a binary code of whether or not it was present (\code{PA}).
  
  \code{covariates}: a \code{RasterBrick} with three layers, each giving the value of a covariate to be used in predicting the distribution of the hypothetical disease across a synthetic landscape. The first two covariates (code{cov_a} and \code{cov_b}) are continuous and the third (\code{cov_c}) is discrete.
  
  \code{consensus}: a \code{RasterLayer} giving the evidence consensus for the disease at regional level in the synthetic landscape. The evidence consensus values are a score ranging from -100 to 100 descrinbing the level of certainty of presence or absence of the disease at a regional level. A  score of -100 denotes absolute certainty of disease absence and 100 denotes absolute certainty of presence. A score of 0 denotes complete uncertainty as to whether the disease is present or absent.
  
  \code{template}: a \code{RasterLayer} acting as a template and land-sea mask for the synthetic landscape.
  
  \code{admin}: a \code{RasterBrick} with 4 layers giving regional codes for admin levels 0 to 3 for the synthetic landscape.
}

\source{
These objects were all generated using the script \code{gen_fake_data.R} which can be downloaded at \url{https://github.com/SEEG-Oxford/seegSDM/tree/master/fake_data}.
}

