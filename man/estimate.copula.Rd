\name{estimate.copula}
\alias{estimate.copula}

\title{Estimation of the structure of HAC
}
\description{
This function estimates the parameter and determines the structure of Hierarchical Archimedean Copulae.
}
\usage{
estimate.copula(X, type = HAC_GUMBEL, method = TAU, epsilon = 0,
agg.method = "mean", margins = NULL, na.rm = FALSE, max.min = TRUE)
}
\arguments{
  \item{X}{
a matrix, which contains the data. The columns should be labeled. The default names \code{X1, ..., Xdim[X]2} are given otherwise.
}
  \item{type}{
  defines the model. It can be selected between \code{HAC_GUMBEL}, \code{AC_GUMBEL}, \code{HAC_CLAYTON}, \code{AC_CLAYTON} and \code{GAUSS}.
}
  \item{method}{
the estimation method. Either ``Kendall's \eqn{\tau}'' denoted by \code{TAU} or ``Quasi Maximum Likelihood'' denoted by \code{ML}.
}
  \item{epsilon}{
  scalar \eqn{\geq 0}. The variables of consecutive nodes are aggregated, if the difference of their dependency parameters is smaller than \code{epsilon}. For a detailed explanation see also \code{\link{aggregate.hac}}.
}
  \item{agg.method}{
if \eqn{\epsilon > 0}, the new dependency parameter can be determined by \code{"mean"}, \code{"min"} or \code{"max"} of the of the two previous; see also \code{\link{aggregate.hac}}.
}
  \item{margins}{
scalar or vector specifying the marginal distributions. They can be determined nonparametrically denoted by \code{"edf"} or in a parametric way , e.g. \code{"norm"}. According to the latter approach the parameter of the distributions are estimated by ML. Basing on these estimates the values of the univariate margins are computed. If the argument is defined as scalar, all margins are computed according to the specification. Otherwise, different margins may be defined, e.g. \code{c("norm", "t", "edf")} for a 3 dimensional sample. Almost all continuous functions of \code{\link{Distributions}} are available. Inappropriate usage of this argument might lead to misspecified margins, e.g. application of \code{"exp"} though the sample contains negative values. If the sample contains the margins yet, the default value is to use.
}
  \item{na.rm}{
boolean. If \code{na.rm = TRUE} and a row of \code{X} contains \code{NA}, the row is removed and not used for the computation. Corresponding warnings are shown.
}
  \item{max.min}{
boolean. If \code{max.min = TRUE} and an element of \code{X} is \eqn{\geq 1} or \eqn{\leq 0}, it is set to \eqn{1-10^{-6}} and \eqn{10^{-6}} respectively.
}
}
\details{
Because negative dependencies cannot be estimated basing on the implemented HAC-types, the correlations are set to \eqn{0}, if variables of X are negatively correlated and Kendall's \eqn{\tau} is chosen as estimation method.
}
\value{
A \env{hac} object is returned. In the case of \code{type=GAUSS} a matrix is returned.
}
\references{
Okhrin, O., Okhrin, Y. and Schmid, W., Determining the structure and estimation of hierarchical Archimedean copulas, \emph{Journal of Econometrics}, under revision.
	
Joe, H. 2005, Asymptotic efficiency of the two-stage estimation method for copula-based models, \emph{Journal of Multivariate Analysis, Volume 94 Issue 2}, 401-419.
	
Genest, C., Ghoudi, K., and Rivest, L. P. 1995, A Semiparametric Estimation Procedure of Dependence Parameters in Multivariate Families of Distributions, \emph{Biometrika} 82, 543-552.
}
\author{
Ostap Okhrin
}

\examples{
M = matrix(c(2, 0, 0, 0, 3, 4, 0, 0, ~X1 + X5, ~X2, ~X3, ~X4),
ncol = 4, byrow = TRUE)
g.model = hac(HAC_GUMBEL, M)

# sample from copula g.model
x = rHAC(1000, g.model)

# in the following case the true model is binary approximated
est.obj = estimate.copula(x, HAC_GUMBEL, TAU, epsilon = 0)
plot(est.obj, index = TRUE)

# consider also the aggregation of the variables
est.obj = estimate.copula(x, HAC_GUMBEL, TAU, epsilon = 0.1)
plot(est.obj, index = TRUE)

# if margins != NULL
data = rmnorm(1000, varcov = matrix(c(5, 2, 0, 2, 20, 0, 0, 0, 40), 
ncol = 3))

result1 = estimate.copula(data, method = ML, margins = "norm")
result2 = estimate.copula(data, method = ML, margins = c("norm",
"norm", "norm"))
}