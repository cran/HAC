\name{phi, phi.inv}
\alias{phi}
\alias{phi.inv}

\title{Generator function
}
\description{The Archimedean generator function and its inverse.
}
\usage{
phi(x, theta, type = HAC_GUMBEL)
phi.inv(x, theta, type = HAC_GUMBEL)
}
\arguments{
  \item{x}{
a scalar, vector or matrix at which the function is evaluated. The support of the functions has to be taken into account, i.e. \eqn{x \in [0, \infty]} for the generator function and \eqn{x \in [0, 1]} for its inverse.
}
  \item{theta}{
  the feasible dependency parameter of the copula, i.e. \eqn{\theta \in [1, \infty)} for the Gumbel and Joe family, \eqn{\theta \in (0, \infty)} for the Clayton and Frank family and \eqn{\theta \in [0, 1)} for the Ali-Mikhail-Haq family.
 
}
  \item{type}{
  the copula model, e.g. \code{\link{HAC_GUMBEL}}.
}
}
\examples{
x = runif(100, min = 0, max = 100)
phi(x, theta = 1.2)

# do not run the following example
# phi(x, theta = 0.8) # In phi(x, theta = 0.8) : theta >= 1  is required.
}
