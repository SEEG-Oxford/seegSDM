\name{accuracy statistics}
\alias{rmse}
\alias{devBern}

\title{
Functions to compute accuracy statistics for SDMs
}

\description{
Given a validation set, calculate root mean-squared error (\code{rmse}) or predictive deviance for Bernoulli (1 or 0) data (\code{devBern}).
}

\usage{
rmse(truth, prediction)
devBern(truth, prediction)
}

\arguments{
  \item{truth}{
  The true value from the validation set.
  \code{rmse} is a general purpose error metric and can theoretically be used for any type of data.
  \code{devBern} is valid only for binary (1, 0) data.
}
  \item{prediction}{
  The predicted value from the model.
  For \code{devBern} this should be the predicted probability of presence (1 in the \code{truth}).
}
}

\value{
A single number giving the metric of interest.
}

\seealso{
\code{\link{dbinom}}
}

\examples{
# generate 100 fake datapoints
n <- 100

# 'true' probability of presence
prob <- runif(n)
# true presence/absence
y <- rbinom(n, 1, prob)

# predicted probability of presence
pred <- runif(n)

# Bernoulli deviance of the prediction from the validation set
devBern(y, pred)

# root mean squared error of predicted from true probability of presence
rmse(prob, pred)
}
