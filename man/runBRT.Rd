\name{runBRT}
\alias{runBRT}

\title{
Run a boosted regression tree model using Sam's default settings
}
\description{
A wrapper to run a BRT model using \code{\link{gbm.step}} (or optionally just \code{\link{gbm}}) with Sam's default parameter settings. Covariate effect curves, relative influences and a prediction map on the probability scale are returned. A function to define regression weights can be specified through \code{wt.fun}.
BRT models sometimes fail to converge and the \code{gbm.step} implementation fails silently, returning \code{NULL}. \code{runBRT} instead attempts to run a model \code{max_tries} times and fails with an error if it still hasn't converged.
To run a BRT model \emph{without} optimising the number of trees you can set \code{step} to \code{FALSE} - which should be much faster.
}

\usage{
runBRT(data, gbm.x, gbm.y,
      pred.raster, gbm.coords = NULL,
      wt.fun = function(PA) rep(1, length(PA)),
      max_tries = 5, verbose = FALSE, tree.complexity = 4, learning.rate = 0.005,
      bag.fraction = 0.75, n.trees = 10, n.folds = 10,
      max.trees = 40000, step.size = 10, step = TRUE, ...)
}

\arguments{
  \item{data}{
    Input dataframe.
  }
  \item{gbm.x}{
    Index for columns containing covariate values.
  }
  \item{gbm.y}{
    Index for column containing presence/absence code (1s or 0s).
  }
  \item{pred.raster}{
    A \code{RasterBrick} or \code{RasterStack} object to predict the model to.
  }
  \item{gbm.coords}{
    Optional index for two columns containing coordinates of records. This is required if you later want to calculate validation statistics using pair-wise distance sampling (setting \code{pwd = TRUE} in \code{\link{getStats}}). Set to \code{NULL} (the default) if not required.
  }
  \item{wt.fun}{
A function to specify the weighting applied to resence or background points. The function must take a vector of 1s and 0s as input and return a vector of the same length giving regression weights. The default function implies no weighting of points. To apply a 50:50 weighting of presence and absence records (mimicking a  prevalence of 0.5) use: \code{wt.fun = function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA))}.
}
  \item{max_tries}{
How many time to try and get gbm.step to converge before throwing an error.
}
  \item{verbose}{
Passed to \code{gbm.step}, whether to report on progress.
}
  \item{tree.complexity}{
Passed to \code{gbm.step}, number of bifurcations in each individual tree.
}
  \item{learning.rate}{
Passed to \code{gbm.step}, how small to shrink the contribution of each tree in the final model
}
  \item{bag.fraction}{
Passed to \code{gbm.step}, proportion of datapoints used in selecting variables
}
  \item{n.trees}{
Passed to \code{gbm.step}, \emph{initial} number of trees to fit. \code{gbm.step} optimises this parameter.
}
  \item{n.folds}{
Passed to \code{gbm.step}, number of folds in each round of cross validation.
}
  \item{max.trees}{
Passed to \code{gbm.step}, maximum number of trees to fit before stopping the stepping algorithm.
}
  \item{step.size}{
Passed to \code{gbm.step}, number of trees to add at each iteration.
}
  \item{step}{
Whether to run the \code{gbm.step} procedure. If \code{FALSE} a single model is run with {n.trees} trees. Note that using the default value of \code{n.trees = 10} is likely to be a bad move.
  }
  \item{\dots}{
Additional functions to pass to \code{gbm.step}.
}
}

\value{
A list containing four elements
  \item{model }{the fitted gbm model}
  \item{effects }{a list of effect curves with one element ofr each covariate}
  \item{relinf }{a vector of relative influence estimates for each covariate}
  \item{pred }{a \code{RasterLayer} giving predictions on the probability scale.}
  \item{coords }{a dataframe giving the coordinates of the training points (or \code{NULL} if \code{gbm.coords = NULL})}
}

\seealso{
\code{\link{gbm.step}}, \code{\link{getRelInf}}, \code{\link{getEffectPlots}}, \code{\link{combinePreds}}
}
\examples{
# TO DO
}