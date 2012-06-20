\name{hac}
\alias{hac}
\alias{hac.full}
\alias{print.hac}

\title{Construction of \env{hac} objects
}
\description{
\env{hac} objects are required as input argument for several functions, e.g. \code{\link{plot.hac}} and \code{\link{rHAC}}. These objects can be constructed by \code{hac} and \code{hac.full}. The latter function produces only fully nested Archimedean copulae, whereas \code{hac} can construct any dependence structure.}
\usage{
hac(type = HAC_GUMBEL, tree = NULL, corr = NULL)
hac.full(type = HAC_GUMBEL, y, theta)
\method{print}{hac}(x, digits = 2, \dots)
}
\arguments{
  \item{y}{
a vector containing the variables, which are denoted by a \code{\link{character}}, e.g. \code{"X1"}.
}
  \item{theta}{
a vector containing the dependency parameters of the HAC. The parameters should be ordered from top to down. The length of \code{theta} must be equal to \code{length(y) - 1}.
}
\item{tree}{
a \code{\link{list}} object of the general structure \code{list(..., numeric(1))}. The last argument of the list, \code{numeric(1)}, denotes the dependency parameter. The arguments \code{...} are either of the same structure or of the class \code{\link{character}}. The character objects denote variables and embedded lists refer to nested subcopulae.
}
\item{corr}{
if \code{type = GAUSS}, \code{corr} should be the correlation matrix.
}
\item{type}{
all copula-types are admissible, e.g. \code{\link{HAC_GUMBEL}}.
}
  \item{x}{
a \env{hac} object.
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
\item{tree}{the structure of the HAC. If \code{type = GAUSS}, the correlation matrix.}
}
\examples{
# it might be helpful to plot the hac objects
# Example 1: 4-dim AC
tree = list("X1", "X2", "X3", "X4", 2)
AC = hac(tree = tree)

# Example 2: 4-dim HAC
y = c("X1", "X4", "X3", "X2")
theta = c(2, 3, 4)

HAC1 = hac.full(y = y, theta = theta)
HAC2 = hac(tree = list(list(list("X2", "X3", 4), 
"X4", 3), "X1", 2))
tree2str(HAC1) == tree2str(HAC2) # [1] TRUE

# Example 3: 9-dim HAC

HAC = hac(tree = list(c("X6", "X5"), list("X2", "X4", "X3", 4.4),
list("X1", "X7", 3.3), list("X8", "X9", 4), 2.3))
plot(HAC)

# Example 4: Gaussian-copula
corr = diag(1, 3)
corr[lower.tri(corr)] = corr[upper.tri(corr)] = seq(0.1, 0.3, by = 0.1)
gauss = hac(type = GAUSS, corr = corr) 
}
