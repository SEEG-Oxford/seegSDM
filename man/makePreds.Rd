\name{makePreds}
\alias{makePreds}

\title{
Using the model output from \code{\link{runBRT}}, makes a prediction to a set of covariates.
}
\description{
Given the output from \code{\link{runBRT}} and a rasterbrick of covariates to predict to, the function makes (and returns) a prediction to the rasters based on the fitted model from \code{\link{runBRT}}. The function will throw an error if the column names for the covariates used in \code{\link{runBRT}} do not match the prediction covariate names.
}
\usage{
makePreds(object,
          pred_covs)
}

\arguments{
  \item{object}{
  Output from \code{\link{runBRT}}.
  }
  \item{pred_covs}{
  A rasterbrick of covariates to predict to.
  }
}

\value{
A raster containing prediction values.
}

\seealso{
\code{\link{runBRT}}
}
\examples{
# fit BRT model
example(runBRT)

# reload the covariate rasters
data(covariates)

# change names so function produces an error 
prediction_covs <- covariates

names(prediction_covs)[names(prediction_covs)=='cov_c'] <- 'cov_d'

# get predictions
\dontrun{
# this should produce an error - the columns are wrongly named
pred1 <- makePreds(model, prediction_covs)
}
pred2 <- makePreds(model, covariates)
}

