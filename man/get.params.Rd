\name{get.params}
\alias{get.params}

\title{Dependency parameters of a HAC
}
\description{
This function returns the dependency parameters. They are ordered from top to down and left to right.
}
\usage{
get.params(hac, sort.v = FALSE, \dots)
}
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
  \item{sort.v}{
boolean. If \code{sort.v = TRUE}, the output is sorted.
}
  \item{\dots}{
further arguments passed to \code{\link{sort}}.
}
}
\seealso{
\code{\link{tree2str}}
}
\examples{
# construct a copula model
tree = list(list("X1", "X5", "X2", 4), list("X3", "X4", "X6", 3), 2)
model = hac(HAC_GUMBEL, tree)

# return the parameter
get.params(model) # [1] 2 4 3
get.params(model, sort.v = TRUE, decreasing = TRUE) # [1] 4 3 2
}