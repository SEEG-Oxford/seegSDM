\name{sdWeighted}
\alias{sdWeighted}

\title{
Calculate weighted standard deviation of a vector.
}
\description{
Given a numeric vector (\code{x}), the corresponding vector of weights (\code{weights}) and the weighted mean of \code{x} (\code{weighted_mean}), calculate the weighted standard deviation of \code{x}.
}
\usage{
sdWeighted(x, weights, weighted_mean)
}

\arguments{
  \item{x}{
The numeric vector for which the weighted standard deviation is to be calculated.
}
  \item{weights}{
  A vector of the same length as \code{x} giving the weights.
}
  \item{weighted_mean}{
  The weighted mean of \code{x}, perhaps calculated using \code{\link{weighted.mean}}.
}
}

\value{
  A numeric vector of length 1 giving the weightyed mean of \code{x}.
}

\seealso{
\code{\link{weighted.mean}}
}
\examples{
x <- rnorm(10)
wt <- runif(10)
mn <- weighted.mean(x, wt)
std <- sdWeighted(x, wt, mn)
}
