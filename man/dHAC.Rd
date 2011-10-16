\name{dAC, dHAC}
\alias{dHAC}
\alias{dAC}
\title{Probability density function
}
\description{
\code{dHAC} returns the values of the pdf for 3-dimensional HAC, whereas the function \code{dAC} returns the values of a pdf for 2-dim AC.
}
\usage{
dHAC(X, hac, margins = NULL, na.rm = FALSE, max.min = TRUE)
dAC(x, y, theta = 1.0, type = AC_GUMBEL, na.rm = FALSE, max.min = TRUE)
}
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
  \item{theta}{
  the dependency parameter of the copula.
}
  \item{type}{
  choose between \code{AC_GUMBEL}, \code{AC_CLAYTON} and \code{GAUSS}.
}
  \item{X}{
a data matrix. The number of columns (\eqn{2} or \eqn{3}) has to coincide with the dimension of the corresponding copula model. \code{X} has to contain at least to rows (observations), because the values of the pdf cannot be computed else. The column names have to equal the names of the variables of the argument \code{hac}.
}
  \item{x, y}{
data vectors.
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
\details{
To do: extending the dimension of \code{dHAC}.
}
\value{
A vector containing the values of the pdf.
}
\references{
Savu, C. and Trede, M. 2010, Hierarchies of Archimedean copulas, \emph{Quantitative Finance} 10, 295-304.
}
\author{
Ostap Okhrin
}
\seealso{
\code{\link{pHAC}}
}
\examples{
# AC example
# the underlying model
g.model = hac(type = AC_GUMBEL, X = 1.5, dim = 2)

# sample from copula g.model
sample = rHAC(100, g.model)

# returns the pdf at each point of the sample
values_1 = dHAC(sample, g.model)
values_2 = dAC(sample[, 1], sample[, 2], theta = 1.5, type = AC_GUMBEL)
round(values_1) == round(values_2) # TRUE


# HAC example
# the underlying model
y = c(~X1, ~X2, ~X3)
theta = c(1.5, 3)
g.model = hac.full(type = HAC_GUMBEL, y, theta)

# sample from copula g.model
sample = rHAC(100, g.model)

# returns the pdf at each point of the sample
values = dHAC(sample, g.model)
}