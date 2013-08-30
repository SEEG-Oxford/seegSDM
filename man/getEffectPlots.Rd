\name{getEffectPlots}
\alias{getEffectPlots}

\title{
Combine effect plots for ensemble BRT runs
}

\description{
Given \code{models}, a list of BRT model bootstraps (each element an output from \code{\link{runBRT}}), return a list of matrices for, and optionally plot (\code{\link{plot}}) the mean effect curve for each covariate in the model with confidence regions bounded by \code{quantiles}. The \dots argument is passed to \code{plot} and allows some customisation of the plotting outputs.  
}

\usage{
getEffectPlots(models, plot = FALSE, quantiles = c(0.025, 0.975), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
A list of BRT model bootstraps, each element being an output from \code{runBRT}.
}
  \item{plot}{
Whether to plot the overall effect curves.
}
  \item{quantiles}{
Quantiles from which to calculate the uncertainty regions.
}
  \item{\dots}{
Additional arguments to be passe to \code{\link{plot}}.
}
}

\value{
A list of matrices, one for each covariate, giving the mean and quantiles of the effect curve. Optionally a plot as a side-effect.
}

\seealso{
\code{\link{response}}, \code{\link{plot}}, \code{\link{runBRT}}, \code{\link{quantile}}
}
\examples{
# TO DO
}