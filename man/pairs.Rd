\name{par.pairs}
\alias{par.pairs}

\title{Parameter of the HAC
}
\description{
This function returns a matrix of dependency parameters. They are pairwise ordered, so that the parameters correspond to the lowest node, at which the variables are joined.
}
\usage{
par.pairs(hac, FUN = NULL, \dots)
}
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
  \item{FUN}{
the parameters of the HAC are returned by default. If \code{FUN = "TAU"}, \code{\link{theta2tau}} is applied to the parameters. \code{FUN} can also be a self-defined \code{\link{function}}.
}
  \item{\dots}{
further arguments passed to \code{FUN}.
}
}
\seealso{
\code{\link{get.params}}
}
\examples{
# construct a copula model
tree = list(list("X1", "X5", "X2", 4), list("X3", "X4", "X6", 3), 2)
model = hac(HAC_GUMBEL, tree)

# returns the pairwise parameter
par.pairs(model)

# Kendall's TAU
par.pairs(model, FUN = "TAU")

# sqrt of the parameter
par.pairs(model, function(r)sqrt(r))
}
