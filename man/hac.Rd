\name{hac}
\alias{hac}
\alias{hac.full}
\alias{print.hac}
\alias{hac2nacopula}
\alias{nacopula2hac}

\title{Construction of \env{hac} objects
}
\description{
\env{hac} objects are required as input argument for several functions, e.g. \code{\link{plot.hac}} and \code{\link{rHAC}}. They can be constructed by \code{hac} and \code{hac.full}. The latter function produces only fully nested Archimedean copulae, whereas \code{hac} can construct arbitrary dependence structures for a given family. Moreover, the functions \code{hac2nacopula} and \code{nacopula2hac} ensure the compatability with the \env{copula} package.}
\usage{
hac(type, tree)
hac.full(type, y, theta)
\method{print}{hac}(x, digits = 2, \dots)
hac2nacopula(x)
nacopula2hac(outer_nacopula)
}
\arguments{
  \item{y}{
a vector containing the variables, which are denoted by a \code{\link{character}}, e.g. \code{"X1"}.
}
  \item{theta}{
a vector containing the HAC parameters, which should be ordered from top to down. The length of \code{theta} must be equal to \code{length(y) - 1}.
}
\item{tree}{
a \code{\link{list}} object of the general structure \code{list(..., numeric(1))}. The last argument of the list, \code{numeric(1)}, denotes the dependency parameter. The arguments \code{...} are either of the same structure or of the class \code{\link{character}}. The character objects denote variables and embedded lists refer to nested subcopulae.
}
\item{type}{
all copula-types are admissible, see \code{\link{phi}} for an overview of implemented families.
}
  \item{x}{
a \env{hac} object.
}
 \item{outer_nacopula}{
an \code{nacopula} object. The variables of the \code{outer_nacopula} object \code{1, 2, ...} are translated into the \code{\link{character}}s \code{"X1", "X2", ...}.
}
  \item{digits}{
specifies the digits, see \code{\link{tree2str}}.
}
  \item{\dots}{
arguments to be passed to \code{\link{cat}}.
}
}
\value{
A \env{hac} object is returned.
\item{type}{the specified copula type.}
\item{tree}{the structure of the HAC.}
}
\references{
Hofert, M. and Maechler, M. 2011, Nested Archimedean Copulas Meet R: The \verb{nacopula} Package, \emph{Journal of Statistical Software}, 39(9), 1-20, \doi{10.18637/jss.v039.i09}.

Hofert, M., Kojadinovic, I., Maechler, M. and Yan, J. 2015, \verb{copula}: Multivariate Dependence with Copulas, \emph{R package version 0.999-14}, \url{https://CRAN.R-project.org/package=copula}.

Kojadinovic, I., Yan, J. 2010, Modeling Multivariate Distributions with Continuous Margins Using the \verb{copula} R Package, \emph{Journal of Statistical Software}, 34(9), 1-20. \doi{10.18637/jss.v034.i09}.

Okhrin, O. and Ristig, A. 2014, Hierarchical Archimedean Copulae: The \verb{HAC} Package", \emph{Journal of Statistical Software}, 58(4), 1-20, \doi{10.18637/jss.v058.i04}.

Yan, J. 2007, Enjoy the Joy of Copulas: With a Package \verb{copula}, \emph{Journal of Statistical Software}, 21(4), 1-21, \doi{10.18637/jss.v021.i04}.
}
\examples{
# it might be helpful to plot the hac objects
# Example 1: 4-dim AC
tree = list("X1", "X2", "X3", "X4", 2)
AC = hac(type = 1, tree = tree)

# Example 2: 4-dim HAC
y = c("X1", "X4", "X3", "X2")
theta = c(2, 3, 4)

HAC1 = hac.full(type = 1, y = y, theta = theta)
HAC2 = hac(type = 1, tree = list(list(list("X2", "X3", 4), 
"X4", 3), "X1", 2))
tree2str(HAC1) == tree2str(HAC2) # [1] TRUE

# Example 3: 9-dim HAC

HAC = hac(type = 1, tree = list("X6", "X5", list("X2", "X4", "X3", 4.4),
list("X1", "X7", 3.3), list("X8", "X9", 4), 2.3))
plot(HAC)
}
