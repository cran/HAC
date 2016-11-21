\name{emp.copula} \alias{emp.copula} \alias{emp.copula.self} \title{Empirical copula } \description{\code{emp.copula} and \code{emp.copula.self} compute the empirical copula for a given sample. The difference between these 
functions is, that \code{emp.copula.self} does not require a matrix \code{u}, at which the function is evaluated.} 

\usage{
emp.copula(u, x, proc = "M", sort = "none", margins = NULL, 
na.rm = FALSE, \dots)
emp.copula.self(x, proc = "M", sort = "none", margins = NULL, 
na.rm = FALSE, \dots)
}

\arguments{
  \item{u}{
a matrix, at which the function is evaluated. According to the dimension of the data matrix \code{x}, it can be a scalar, a vector or a matrix. The entries of \code{u} should be within the interval \eqn{[0, 1]}. }  \item{x}{ 
denotes the matrix of marginal distributions, if \code{margins = NULL}. The number of columns should be equal the dimension \eqn{d}, whereas the number of rows should be equal to the number of observations 
\eqn{n}, with \eqn{n > d}. } 
  \item{proc}{
enables the user to choose between two different methods. It is recommended to use the default method, \code{"M"}, because it takes only a small fraction of the computational time of method \code{"A"}. However, method 
\code{"M"} is sensitive with respect to the size of the working memory and therefore, non-applicable for very large datasets.} 
  \item{sort}{
defines, whether the output is ordered. \code{sort = "asc"} refers to ascending values, which might be interesting for plotting and \code{sort = "desc"} refers to descending values.
}
  \item{margins}{
  specifies the margins. The data matrix is assumed to contain the values of the marginal distributions by default, i.e. \code{margins = NULL}. If raw data are used, the margins can be determined nonparametrically, \code{"edf"}, or in parametric way, e.g. \code{"norm"}. See \code{\link{estimate.copula}} for a detailed explanation.
}
  \item{na.rm}{
boolean. If \code{na.rm = TRUE}, missing values, \code{NA}, contained in \code{x} and \code{u} are removed.
}
  \item{\dots}{
arguments to be passed to \code{\link{na.omit}}.
}
}
\details{
The estimated copula follows the formula \deqn{\widehat{C} \left(u_{1}, \dots, u_{d} \right) = n^{-1} \sum_{i=1}^{n} \prod_{j=1}^{d} \mathbf{I} \left\{ \widehat{F}_{j} \left( X_{ij} \right) \leq u_{j} \right\},} where \eqn{\widehat{F}_{j}} denotes the empirical marginal distribution function of variable \eqn{X_{j}}.
}
\value{
A vector containing the values of the empirical copula.
}
\references{
Okhrin, O. and Ristig, A. 2014, Hierarchical Archimedean Copulae: The \verb{HAC} Package", \emph{Journal of Statistical Software}, 58(4), 1-20, \url{https://www.jstatsoft.org/v58/i04/}.
}
\seealso{
\code{\link{pHAC}}
}
\examples{
v = seq(-4, 4, 0.05)
X = cbind(matrix(pt(v, 1), 161, 1), matrix(pnorm(v), 161, 1))

# both methods lead to the same result 
z = emp.copula.self(X, proc = "M") 
which(((emp.copula.self(X[1:100, ], proc = "M") - emp.copula.self(X[1:100, ],
proc = "A")) == 0) == "FALSE")
# integer(0)

# the contour plot
out = outer(z, z)
contour(x = X[,1], y = X[,2], out, main = "Contour Plot", 
xlab = "Cauchy Margin", ylab = "Standard Normal Margin", 
labcex = 1, lwd = 1.5, nlevels = 15)
}
