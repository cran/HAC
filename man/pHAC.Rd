\name{pHAC}
\alias{pHAC}
\alias{pAC}

\title{Cumulative distribution function
}
\description{
This function computes the cdf for a given sample and copula model.
}
\usage{
pHAC(X, hac, margins = NULL, na.rm = FALSE, max.min = TRUE)
}
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
  \item{X}{
a data matrix. The sample has to contain at least \eqn{2} rows (observations), since the values of the cdf cannot be computed otherwise. The column names have to be identical to the names of the variables of the \code{hac} object.
}
  \item{margins}{
scalar or vector specifying how the margins are to compute. They can be determined nonparametrically denoted by \code{"edf"} or in parametric way , e.g. \code{"norm"}. See \code{\link{estimate.copula}} for a detailed explanation.
}
  \item{na.rm}{
boolean. If \code{na.rm = TRUE} and a row of \code{X} contains \code{NA}, the row is removed and not used for the computation. Corresponding warnings are shown.
}
  \item{max.min}{
boolean. If \code{max.min = TRUE} and an element of \code{X} is \eqn{\geq 1} or \eqn{\leq 0}, it is set to \eqn{1-10^{-6}} and \eqn{10^{-6}} respectively.
}
}
\value{
A vector containing the values of the cdf.
}
\references{
Joe, H. 1997, Multivariate Models and Dependence Concepts, \emph{Chapman & Hall}.
	
Nelsen, R.,B. 2006, An Introduction to Copulas, \emph{Spinger}, 2nd Edition.
}
\author{
Ostap Okhrin
}
\seealso{
\code{\link{dHAC}}, \code{\link{dAC}}
}
\examples{
# a hac obejct is defined
M = matrix(c(2, 0, 0, 0, 3, 4, 0, 0, ~X1 + X5, ~X2, ~X3, ~X4 + X6),
ncol = 4, byrow = TRUE)
g.model = hac(HAC_GUMBEL, M)

# sample from copula g.model
sample = rHAC(100, g.model)

# cdf at each vector of the sample
values = pHAC(sample, g.model)
}