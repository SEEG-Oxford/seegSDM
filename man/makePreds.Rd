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
  Output from runBRT.
  }
  \item{pred_covs}{
  A rasterbrick of covariates to predict to.
  }
}

}
\value{
A raster containing prediction values.
}

\seealso{
\code{\link{runBRT}}
}
\examples{
# load the data
data(occurrence)

# load the covariate rasters
data(covariates)

prediction_covs <- covariates

names(prediction_covs)[names(prediction_covs)=='cov_c'] <- 'cov_d'

# load evidence consensus layer
data(consensus)

background <- bgSample(consensus,
                       n= 100,
                       replace=FALSE,
                       spatial=FALSE)

colnames(background) <- c('Longitude', 'Latitude')
background <- data.frame(background)

# combine the occurrence and background records
dat <- rbind(cbind(PA = rep(1, nrow(occurrence)),
                   occurrence[, c('Longitude', 'Latitude')]),
             cbind(PA = rep(0, nrow(background)),
                   background[ ,c('Longitude', 'Latitude')]))

# extract covariate values for each data point
dat_covs <- extract(covariates, dat[, c('Longitude', 'Latitude')])

# combine covariates with the other info
dat_all <- cbind(dat, dat_covs)

model <- runBRT(dat_all,
                gbm.x = 4:6,
                gbm.y = 1,
                n.folds = 5)

# get predictions
pred1 <- makePreds(model, prediction_covs)
pred2 <- makePreds(model, covariates)
}

