\name{getRelInf}
\alias{getRelInf}

\title{
Summarize relative influence statistics for ensemble BRT runs
}
\description{
Given \code{models}, a list of BRT model bootstraps (each element an output from \code{\link{runBRT}}), return a matrix of, and optionally boxplot (\code{\link{boxplot}}) the mean and quantiles (determined by \code{quantiles}) of the relative influence of each covariate in the BRT ensembles. The \dots argument is passed to \code{boxplot} and allows some customisation of the plotting outputs.
}
\usage{
getRelInf(models, plot = FALSE, quantiles = c(0.025, 0.975), ...)
}

\arguments{
  \item{models}{
A list of BRT model bootstraps, each element being an output from \code{runBRT}.
}
  \item{plot}{
Whether to plot a box and whisker plot.
}
  \item{quantiles}{
Quantiles from which to calculate the uncertainty regions.
}
  \item{\dots}{
Additional arguments to be passe to \code{boxplot}.
}
}

\value{
A matrix giving the means and quantiles of relative influence statistics for each covariate int he BRT ensemble.
}

\seealso{
\code{\link{relative.influence}}, \code{\link{runBRT}}, \code{\link{boxplot}}
}
\examples{
# TO DO
}
