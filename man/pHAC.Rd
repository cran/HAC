\name{pHAC}
\alias{pHAC}
\alias{pAC}

\title{Cumulative distribution function
}
\description{
This function computes the cdf for a given sample and copula model.
}
\usage{
pHAC(X, hac)
}
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
  \item{X}{
a data matrix. The sample has to contain at least \eqn{2} rows (observations), since the values of the cdf cannot be computed otherwise. The column names have to be identical to the names of the variables of the \code{hac} object.
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