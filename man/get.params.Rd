\name{get.params}
\alias{get.params}

\title{Dependency parameters of a HAC
}
\description{
This function returns the dependency parameters. They are ordered from top to down and left to right.
}
\usage{
get.params(hac)
}
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
}
\author{
Ostap Okhrin
}
\seealso{
\code{\link{tree2str}}
}
\examples{
# construct a copula model
M = matrix(c(2, 0, 0, 0, 4, 3, 0, 0, ~X1 + X5, ~X2, ~X3, ~X4 + X6),
ncol = 4, byrow = TRUE)
g.model = hac(HAC_GUMBEL, M)

# return the parameter
get.params(g.model) # [1] 2 4 3
}