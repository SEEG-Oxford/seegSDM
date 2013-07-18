\name{seegSDM-package}
\alias{seegSDM-package}
\alias{seegSDM}
\docType{package}
\title{
Streamlined functions for species distribution modelling, as applied in seeg research group.
}
\description{
The present focus is on ensemble boosted regression tree methods such as those used in Bhatt et al. (2013).
This will be developed to incorporate other statistical models, methods for selection of pseudo-absence data and other modelling choices.
}
\details{
\tabular{ll}{
Package: \tab seegSDM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-07-17\cr
License: \tab What license is it under?\cr
}
Functions are provided to assist with the seeg SDM workflow:
  \begin{enumerate}
  \item quality control of occurrence and covariate data
  \item generation of pseudo-absence data
  \item fitting ensembles of models in parallel
  \item combining model ensembles to produce final predictions, uncertainty estimates and model summaries
  \end{enumerate}
}
\author{
Nick Golding, based on functions by Samir Bhatt
}
\references{
Bhatt et al. (2013)  The global distribution and burden of dengue. Nature \url{http://dx.doi.org/10.1038/nature12060}
}