\name{copMult}
\alias{copMult}

\title{2dim copula
}
\description{
This function returns the values for \eqn{d} dimensional Archimedean copulae.
}
\usage{
copMult(X, theta, type)
}
\arguments{
  \item{X}{
a \eqn{n \times d} matrix, where \eqn{d} denotes the amount of variables.
}
  \item{theta}{
the dependency parameter of the copula.
}
  \item{type}{
  defines the model. It can be selected between \code{HAC_GUMBEL}, \code{AC_GUMBEL}, \code{HAC_ROTATED_GUMBEL}, \code{HAC_CLAYTON} and \code{AC_CLAYTON}.
}
}
\details{
If warnings are returned, see \code{\link{phi}}.
}
\value{
A numeric vector containing the values of the copula.
}
\author{
Ostap Okhrin
}
\seealso{
\code{\link{pHAC}}}
\examples{
# the arguments are defined
x = runif(100)
X = cbind(x, x, x)

# the values are computed
cop = copMult(X, theta = 1.5, type = HAC_GUMBEL)
}