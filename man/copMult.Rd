\name{copMult}
\alias{copMult}

\title{\eqn{d}-dim copula
}
\description{
This function returns the values for \eqn{d}-dimensional Archimedean copulae.
}
\usage{
copMult(X, theta, type)
}
\arguments{
  \item{X}{
a \eqn{n \times d} matrix, where \eqn{d} refers to the dimension of the copula.
}
  \item{theta}{
the parameter of the copula.
}
  \item{type}{
all copula-types produced by Archimedean generators, see \code{\link{phi}} for an overview of implemented families.
}
}
\details{
If warnings are returned, see \code{\link{phi}}.
}
\value{
A vector containing the values of the copula.
}
\seealso{
\code{\link{pHAC}}}
\examples{
# the arguments are defined
X = matrix(runif(300), ncol = 3)

# the values are computed
cop = copMult(X, theta = 1.5, type = 1)
}
