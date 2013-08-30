\name{subsample}
\alias{subsample}

\title{
Sample a random subset of a dataframe containing presence and absence data
}
\description{
Extract a random subset of \code{n} records from \code{data}, ensuring that there are at least \code{minimum} presence and absence records in the subset. The presence/absence column is identified by \code{prescol}.
}
\usage{
subsample(data, n, minimum = c(5, 5), prescol = 1, replace = FALSE)
}

\arguments{
  \item{data}{
A dataframe or matrix.
}
  \item{n}{
The number of records required in the subset.
}
  \item{minimum}{
An integer vector of length two giving the minimum number of presence (element one) and absence (element two) records in the subset.
}
  \item{prescol}{
An integer giving the column number for the presence/absence code (containing 1s for presences and 0s for absences.)
}
  \item{replace}{
Whether to sample with replacement.
}
}
\value{
A dataframe or matrix (whichever \code{data} was).
}
\seealso{
\code{\link{sample}}
}

\examples{
n <- 100
dat <- data.frame(y = rbinom(n, 1, 0.2),
                  x = rnorm(n),
                  z = rnorm(n))
sub <- subsample(dat, 10, minimum = c(3, 3), prescol = 1)
sub
}