\name{emp.copula}
\alias{emp.copula}
\alias{emp.copula.self}
\title{Empirical copula
}
\description{
The functions \code{emp.copula} and \code{emp.copula.self} compute the empirical copula for a given sample. The difference between these functions is, that \code{emp.copula.self} takes the sample itself as grid \code{u}, such that the grid has not to be specified manually.
}

\usage{
emp.copula(u, x, proc = "M", sort = "none", na.rm = FALSE, max.min = TRUE)
emp.copula.self(x, proc = "M", sort = "none", na.rm = FALSE, max.min = TRUE)
}

\arguments{
  \item{u}{
the grid for the computation of the empirical copula. It can be a scalar, a vector or a matrix. The number of columns has to be equal to the number of columns of the matrix \code{x}. The values of \code{u} should be within the interval \eqn{[0, 1]}. 
}  \item{x}{
denotes the matrix of variables, which are marginal distributions in practice. The values should be within the interval \eqn{[0, 1]}. The number of columns should correspond to the number of variables, whereas the number of rows should be equal to the number of observations.
}
  \item{proc}{
enables the user to choose between two different methods. Generally, it is advised to choose \code{"M"} (given by default), because it takes only \eqn{1/25} of the computational time of method \code{"A"}. However, in the case of large datasets it might occur, that method \code{"M"} does not work, since it requires more working memory. When this occurs, the user is advised to set \code{proc = "A"}, since the procedure is less sensitive regarding the size of the working memory.
}
  \item{sort}{
defines, whether the output is ordered. \code{sort = "asc"} refers to ascending values, which might be interesting for plotting and \code{sort = "desc"} refers to descending values.
}
  \item{na.rm}{
boolean. If \code{na.rm = TRUE} and a row of \code{x} or \code{u} contains \code{NA}, the row is removed and not used for the computation. Corresponding warnings are shown.
}
  \item{max.min}{
boolean. If \code{max.min = TRUE} and an element of \code{x} or \code{u} is \eqn{> 1} or \eqn{< 0}, it is set to \eqn{1} and \eqn{0} respectively.
}
}
\details{
The estimated copula follows the formula \deqn{\widehat{C} \left(u_{1}, \dots, u_{d} \right) = n^{-1} \sum_{i=1}^{n} \prod_{j=1}^{d} \mathbf{I} \left\{ \widehat{F}_{j} \left( X_{ij} \right) \leq u_{j} \right\},} where \eqn{\widehat{F}_{j}} denotes the empirical marginal distribution function of variable j. So it is no smoothing implemented yet.
}
\value{
A vector containing the values of the empirical copula.
}
\author{
Ostap Okhrin, Alexander Ristig
}

\seealso{
\code{\link{pHAC}}
}
\examples{
v = seq(-4, 4, 0.005)
X = cbind(matrix(pt(v, 1), 1601, 1), matrix(pnorm(v), 1601, 1))

# both methods lead to the same result 
z = emp.copula.self(X, proc = "M") 
which(((emp.copula.self(X[1:100, ], proc = "M") - emp.copula.self(X[1:100, ],
proc = "A")) == 0) == "FALSE")
# integer(0)

# furthermore recognize the contour plot
out = outer(z, z)
contour(x = X[,1], y = X[,2], out, main = "Contour Plot", 
xlab = "Cauchy Margin", ylab = "Standard Normal Margin", 
labcex = 1, lwd = 1.5, nlevels = 15)
}