\name{getConditionalEffectPlots}
\alias{getConditionalEffectPlots}

\title{
Combine conditional effect plots for ensemble BRT runs
}
\description{
Given a list of BRT model bootstraps (each element an output from \code{\link{runBRT}}), return a list of matrices for, and optionally plot (\code{\link{plot}}), the mean conditional effect curve for each covariate in the model with confidence regions bounded by \code{quantiles}. There is also an option to specify (using \code{hold}) a covariate \code{value} for which conditional effect curves are generated. The \dots argument is passed to \code{plot} and allows some customisation of the plotting outputs. 
}
\usage{
getConditionalEffectPlots(models, plot = FALSE, quantiles = c(0.025, 0.975), hold = NULL, value = NULL, ...)
}

\arguments{
  \item{models}{
A list of BRT model bootstraps, each element being an output from \code{runBRT}.
}
  \item{plot}{
Whether to plot the overall conditional effect curves.
}
  \item{quantiles}{
Quantiles from which to calculate the uncertainty regions.
}
  \item{hold}{
Option to specify the column number of a covariate which is to be held at a particular \code{value} or level (if discrete) to generate the conditional effect curves.
}
  \item{value}{
Value or level (if discrete) at which to hold a specified covariate. 
}
  \item{\dots}{
Additional arguments to be passed to \code{\link{plot}}.
}
}

\value{
A list of matrices, one for each covariate, giving the mean and quantiles of the conditional effect curve as well as the conditional effect curves for each submodel. Optionally a plot as a side-effect.
}

\seealso{
\code{\link{plot}}, \code{\link{runBRT}}, \code{\link{quantile}}, \code{\link{getEffectPlots}}
}
\examples{

# load the data
data(occurrence)

# load the covariate rasters
data(covariates)

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

# let runBRT know that cov_c is a discrete variable 
dat_all$cov_c <- factor(dat_all$cov_c)

# get random bootstraps of the data (minimum 5 pres/5 abs)
data_list <- replicate(4,
                       subsample(dat_all,
                                 nrow(dat_all),
                                 replace = TRUE,
                                 minimum = c(5, 5)),
                       simplify = FALSE)

model_list <- sfLapply(data_list,
                       runBRT,
                       gbm.x = 4:6,
                       gbm.y = 1,
                       n.folds = 5)

effects <- getConditionalEffectPlots(model_list,
                                     plot = TRUE)

effects2 <- getConditionalEffectPlots(model_list,
                                     plot = TRUE, hold = 2, value = 2)

}


