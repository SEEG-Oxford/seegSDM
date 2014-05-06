\name{runABRAID}
\alias{runABRAID}

\title{
Run a disease risk-mapping workflow for a single disease
}
\description{
A wrapper to produce a disease risk map from disease occurrence data and relevant environmental covariates following a standardised workflow.
This function is intended for use in the ABRAID mapping platform, rather than for bespoke risk maps. As such the methods contained within are likely to change regularly and are not detailed here. At present they are based on the BRT ensemble approach in the Bhatt \emph{et al.} (2013) dengue mapping paper.
}

\usage{
runABRAID(occurrence_path, extent_path, admin_path,
          covariate_path, discrete = rep(FALSE, length(covariate_path)),
          verbose = TRUE, max_cpus = 32)
}

\arguments{
  \item{occurrence_path}{
    A character vector of length one specifying the file path of a csv file containing occurrence data. this file must contain the columns: 'Longitude', 'Latitude' (giving the coordinates of points), 'Weight' (giving the degree of weighting to assign to each occurrence record), 'Admin' (giving the admin level of the record - e.g. 1, 2 or 3 for polygons or -999 for points) and 'GAUL' (the GAUL code corresponding to the admin unit for polygons, or NA for points).
  }
  \item{extent_path}{
    A character vector of length one specifying the file path of an  ASCII raster file giving the definitive extents of the disease.
  }
  \item{admin_path}{
    A character vector specifying the file paths of ASCII raster files giving the GAUL codes for administrative regions at levels one and above. These must be given in ascending order starting at 1.
  }
  \item{covariate_path}{
    A character vector specifying the file paths of ASCII raster files giving the environmetnal covariates used to predict risk of the disease.
  }
  \item{discrete}{
    A logical vector of the same length as \code{covariate_path} stating whether each covariate should be considered as continuous (\code{FALSE}) or discrete (\code{TRUE}).
  }
  \item{verbose}{
  Whether to return messages from functions called in the script
  }
  \item{max_cpus}{
  The maximum number of CPUs to use when running the ensemble. At present \code{runABRAID} runs 64 bootstrap submodels per disease, so the number of cpus used in the cluster will be set at \code{min(64, max_cpus)}.
  }
}

\value{
Nothing is returned, but the following files are written to the working directory:
  \item{mean_prediction.asc }{a raster map giving the mean predicted disease risk}
  \item{prediction_uncertainty.asc }{a raster map giving uncertainty in the mean predicted disease risk (specifically the range of the 95\% confidence interval)}
  \item{statistics.csv }{a csv giving the cross-validation statistics for each submodel in the ensemble}
  \item{relative_influence.csv }{a csv file summarizing the relative influence of each covaraite on the final model}
  \item{effect_plots.png }{a multi-panel plot of the marginal effect curves for each of the covariates used in the model}
}
