\name{getStats}
\alias{getStats}

\title{
Get model validation statistics for a model object created by \code{\link{runBRT}}.
}
\description{
Given an object returned by \code{runBRT}, extract \code{\link{Kappa}}, \code{\link{auc}}, \code{\link{sensitivity}} and \code{\link{specificity}} and proportion correctly classified (\code{\link{pcc}}) validation statistics - all calculated using the \code{\link{PresenceAbsence}} package. The estimates returned are the means of the validation statistics and their standard deviations calculated against the witheld data for each of the \code{n.folds} folds in the BRT run. That is, if the arguments \code{n.folds = 10, bag.fraction = 0.75} are passed to \code{runBRT}, the resulting BRT model will be an average of 10 separate BRT models, each of them train on a subset of 75\% of the data. The validation statistics for each fold are calculated by comparing the predictions of each fold model against the 25\% of the data which was witheld for that fold. Estimated standard deviations for these statistics are also calculated (by the functions in the \code{PresenceAbsence} package). The mean of these statistics across the 10 folds is what is reported.

This function does not yet calculate these statistics using the pairwise distance sampling procedure (\code{\link{pwdSample}}) of Hijmans, though future versions will.

}
\usage{
getStats(object)
}

\arguments{
  \item{object}{
A list of BRT model bootstraps, each element being an output from runBRT.
}
}

\value{
A vector giving the mean cross-validation statistics and mean standard deviations for these across the folds (see decription for details).
}

\seealso{
\code{\link{runBRT}}, \code{\link{PresenceAbsence}}, \code{\link{Kappa}}, \code{\link{auc}}, \code{\link{sensitivity}}, \code{\link{specificity}}, \code{\link{pcc}}, \code{\link{pwdSample}}
}
\examples{
# TO DO
}
