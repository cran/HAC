\name{definitions}
\alias{AC_CLAYTON}
\alias{AC_GUMBEL}
\alias{ML}
\alias{TAU}
\alias{HAC_CLAYTON}
\alias{HAC_GUMBEL}
\alias{HAC_ROTATED_GUMBEL}
\alias{GAUSS}
\alias{theta.eps}

\title{Definitions
}
\description{In order that the package works appropriately some definitions are needed, which are loaded with the package.
}
\arguments{
  \item{HAC_GUMBEL = 0}{
  if the user is asked for entering a copula type (or model), it is possible to choose between \code{HAC_GUMBEL} ... \code{GAUSS} or the corresponding numbers.}
  \item{AC_GUMBEL = 1}{
  see \code{HAC_GUMBEL}.}
  \item{HAC_ROTATED_GUMBEL = 2}{
  see \code{HAC_GUMBEL} . Furthermore a Rotated Gumbel copula is defined as Gumbel copula rotated by \eqn{180^\circ}.}
  \item{HAC_CLAYTON = 3}{
  see \code{HAC_GUMBEL}.}
  \item{AC_CLAYTON = 4}{
  see \code{HAC_GUMBEL}.}
  \item{GAUSS = 5}{
  see \code{HAC_GUMBEL}.}
    \item{TAU = 0, ML = 1}{
  there are two procedures regarding the estimation of HAC. It is up to the user to choose between Kendalls's \eqn{\tau} and Quasi Maximum Likelihood.
}
  \item{theta.eps = 0.001}{
  a scalar, which is called during the estimation procedure. The parameter space of the QML-procedure is constraint by the inverse of Kendall's rank correlation coefficient of the interval \eqn{(\theta_{eps}, 1 - \theta_{eps})}. If necessary the parameter space can be adjusted by changing the value of \code{theta.eps}.
}
}