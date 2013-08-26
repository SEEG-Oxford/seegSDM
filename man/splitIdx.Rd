\name{splitIdx}
\alias{splitIdx}

\title{
Get a list of indices to split up a vector for parallel processing
}
\description{
Create a list of indices to split up a vector \code{n} into chunks with a length of at most \code{maxn}.
}
\usage{
splitIdx(n, maxn = 1000)
}

\arguments{
  \item{n}{
  Length of the vector to be split up.
}
  \item{maxn}{
  The maximum size (number of records) of each chunk.
}
}

\value{
A list, each element containing a vector of lenngth two giving the start and end indices for a chunk.
}


\examples{
idx <- splitIdx(n = 26, maxn = 3)

# use lapply on index to plot corresponding letters
lapply(idx, function(idx) letters[idx[1]:idx[2]])
}
