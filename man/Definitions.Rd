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
\description{
Global parameters are defined. If the user is asked for a copula type, choose between \code{HAC_GUMBEL} ... \code{GAUSS} or the corresponding numbers.
}
\arguments{
  \item{HAC_GUMBEL = 0}{--}
  \item{AC_GUMBEL = 1}{--}
  \item{HAC_ROTATED_GUMBEL = 2}{Gumbel copula rotated by \eqn{180^\circ}.}
  \item{HAC_CLAYTON = 3}{--}
  \item{AC_CLAYTON = 4}{--}
  \item{GAUSS = 5}{--}
    \item{TAU = 0, ML = 1}{
  there are two procedures regarding the estimation of HAC. It can be chosen between Kendalls's \eqn{\tau} and Quasi Maximum Likelihood.
} 
}
