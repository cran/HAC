\name{theta2tau, tau2theta}
\alias{theta2tau}
\alias{tau2theta}

\title{Kendall's rank correlation coefficient
}
\description{
Kendall's rank correlation coefficient and its inverse.
}
\usage{
theta2tau(theta, type)
tau2theta(tau, type)
}
\arguments{
  \item{theta}{the dependency parameter. It can be either a scalar, a vector or a matrix and has to lie within a certain interval, i.e. \eqn{\theta \in [1, \infty)} for the Gumbel and Joe family, \eqn{\theta \in (0, \infty)} for the Clayton and Frank family and \eqn{\theta \in [0, 1)} for the Ali-Mikhail-Haq family.
}
\item{tau}{Kendall's rank correlation coefficient. It can be either a scalar, a vector or a matrix and it is to ensure, that \eqn{\tau \in [0,1)}.
}
  \item{type}{all types are available, see \code{\link{phi}} for an overview of implemented families.
}
}
\examples{
# computation of the dependency parameter
x = runif(10)
theta = tau2theta(x, type = 1)

# computation of kendall's tau
y = runif(10, 1, 100)
tau = theta2tau(y, type = 1)
}
