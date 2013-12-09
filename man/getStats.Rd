\name{getStats}
\alias{getStats}

\title{
Get model validation statistics for a model object created by \code{\link{runBRT}}.
}
\description{
Given an object returned by \code{runBRT}, extract \code{\link{devBern}}, \code{\link{rmse}}, \code{\link{auc}}, \code{\link{Kappa}}, \code{\link{sensitivity}} and \code{\link{specificity}} and proportion correctly classified (\code{\link{pcc}}) validation statistics - calculated using either the \code{\link{PresenceAbsence}} or \code{seegSDM} functions. Note that \code{auc} is calculated with a \code{seegSDM} clone of the auc function in \code{PresenceAbsence} but in which worse-than-random AUC scores are not inverted.

If \code{cv = TRUE} the estimates returned are the means of the validation statistics and their standard deviations calculated against the witheld data for each of the \code{n.folds} folds in the BRT run. That is, if the arguments \code{n.folds = 10, bag.fraction = 0.75} are passed to \code{runBRT}, the resulting BRT model will be an average of 10 separate BRT models, each of them train on a subset of 75\% of the data. The validation statistics for each fold are calculated by comparing the predictions of each fold model against the 25\% of the data which was witheld for that fold. Estimated standard deviations for these statistics are also calculated (by the functions in the \code{PresenceAbsence} package). The mean of these statistics across the 10 folds is what is reported.

If \code{cv = TRUE} \emph{and} \code{pwd = TRUE}, these cv statistics are calculated using the pairwise distance sampling procedure (\code{\link{pwdSample}}) of Hijmans (2012) to ensure that accuracy statistics are not inflated by increasing the pseudo-absence selection distance in presence-only models

If \code{cv = FALSE} the statistics are calculated once on the full training set using the final full model. \code{pwd = TRUE} cannot be used in this case.
}

\usage{
getStats(object, cv = TRUE, pwd = TRUE, threshold = 1, ...)
}

\arguments{
  \item{object}{
A list of BRT model bootstraps, each element being an output from runBRT.
}
  \item{cv}{
Whether to calculate cross-validation statistics using folds (\code{cv = TRUE}) or the training-set validation statistics from the final model(\code{cv = FALSE}).
}

  \item{pwd}{
Whether to use the pairwise distance sampling procedure (\code{\link{pwdSample}}) of Hijmans (2012) to ensure that accuracy statistics are not inflated by increasing the pseudo-absence selection distance in presence-only models. Note that this procedure can only be used with \code{cv = TRUE} and where coordinates are available in object (see argument \code{gbm.coords} to \code{runBRT}).
}

  \item{threshold}{
The threshold distance for the pairwise distance sampling procedure, passed directly to the argument \code{tr} in \code{pwdSample}. This will be ignored if \code{pwd = FALSE}.
}

  \item{\dots}{
Other arguments to pass to \code{pwdSample}, such as\code{lonlat}, which you may want to set if the coordinates are not longitude/latitude. Note that the \code{n} argument is fixed at one.
}
}

\value{
A vector giving the mean cross-validation statistics and mean standard deviations for these across the folds (see decription for details).
}

\references{
  Hijmans, R.J., 2012. Cross-validation of species distribution models: removing spatial sorting bias and calibration with a null-model. Ecology 93: 679-688
}

\seealso{
\code{\link{runBRT}}, \code{\link{PresenceAbsence}}, \code{\link{Kappa}}, \code{\link{auc}}, \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{pcc}}, \code{\link{pwdSample}}
}
\examples{
# TO DO
}
