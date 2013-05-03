\name{definitions}
\alias{AC_CLAYTON}
\alias{AC_GUMBEL}
\alias{ML}
\alias{TAU}
\alias{FML}
\alias{RML}
\alias{HAC_CLAYTON}
\alias{HAC_GUMBEL}
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
  \item{HAC_CLAYTON = 3}{--}
  \item{AC_CLAYTON = 4}{--}
  \item{GAUSS = 5}{--}
    \item{TAU = 0, ML = 1, FML = 2, RML = 3}{
  there are four procedures regarding the estimation of HAC. It can be chosen between Kendalls's \eqn{\tau}, Quasi Maximum Likelihood, Full Maximum Likelihood and Recursive Maximum Likelihood.
} 
}
