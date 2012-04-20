\name{estimate.copula}
\alias{estimate.copula}

\title{Estimation of the structure of HAC
}
\description{
This function estimates the parameter and determines the structure of Hierarchical Archimedean Copulae.
}
\usage{
estimate.copula(X, type = HAC_GUMBEL, method = ML, epsilon = 0, agg.method = "mean", 
margins = NULL, theta.eps = 0.001, na.rm = FALSE, max.min = TRUE, \dots)
}
\arguments{
  \item{X}{
a \eqn{n \times d} matrix. The columns should be labeled. The default names \code{X1, X2, X3, ...} are given otherwise.
}
  \item{type}{
  defines the model. Select between \code{HAC_GUMBEL}, \code{AC_GUMBEL}, \code{HAC_CLAYTON} and \code{AC_CLAYTON}.
}
  \item{method}{
the estimation method. Either ``Quasi Maximum Likelihood'' denoted by \code{ML} or ``Kendall's \eqn{\tau}'' denoted by \code{TAU}.
}
  \item{epsilon}{
  scalar \eqn{\geq 0}. The variables of consecutive nodes are aggregated, if the difference of the dependency parameters is smaller than \code{epsilon}. For a detailed explanation see also \code{\link{aggregate.hac}}.
}
  \item{agg.method}{
if \eqn{\epsilon > 0}, the new dependency parameter can be determined by \code{"mean"}, \code{"min"} or \code{"max"} of the two parameters, see \code{\link{aggregate.hac}}.
}
  \item{margins}{
  specifies the margins. The data matrix is assumed to contain the values of the marginal distributions by default, i.e. \code{margins = NULL}. If raw data are used, the margins can be determined nonparametrically, \code{"edf"}, or in a parametric way, e.g. \code{"norm"}. Following the latter approach, the parameters of the distributions are estimated by ML. Basing on these estimates the values of the univariate margins are computed. If the argument is defined as scalar, all margins are computed according to this specification. Otherwise, different margins can be defined, e.g. \code{c("norm", "t", "edf")} for a 3-dimensional sample. Almost all continuous functions of \code{\link{Distributions}} are available. Inappropriate usage of this argument might lead to misspecified margins, e.g. application of \code{"exp"} even though the sample contains negative values.
}
\item{theta.eps}{
a scalar. The parameter space of the estimation procedure is constraint from below by \code{\link{tau2theta}(theta.eps)}.
}
  \item{na.rm}{
boolean. If \code{na.rm = TRUE}, missing values, \code{NA}, contained in \code{X} are removed.
}
  \item{max.min}{
boolean. If \code{max.min = TRUE} and an element of \code{X} is \eqn{\geq 1} or \eqn{\leq 0}, it is set to \eqn{1-10^{-6}} and \eqn{10^{-6}} respectively.
}
  \item{\dots}{
further arguments passed to or from other methods, e.g. \code{\link{na.omit}}.
}
}
\details{
If variables are negatively correlated and Kendall's \eqn{\tau} is chosen as estimation method, the correlations are set to \code{theta.eps}, since the implemented HAC-types does not permit negative dependence. 
}
\value{
A \env{hac} object is returned.
}
\references{
Okhrin, O., Okhrin, Y. and Schmid, W., Determining the structure and estimation of hierarchical Archimedean copulas, \emph{Journal of Econometrics}, under revision.
	
Joe, H. 2005, Asymptotic efficiency of the two-stage estimation method for copula-based models, \emph{Journal of Multivariate Analysis, Volume 94 Issue 2}, 401-419.
	
Genest, C., Ghoudi, K., and Rivest, L. P. 1995, A Semiparametric Estimation Procedure of Dependence Parameters in Multivariate Families of Distributions, \emph{Biometrika} 82, 543-552.
}
\examples{
# define the copula model
tree = list(list("X1", "X5", 3), list("X2", "X3", "X4", 4), 2)
model = hac(tree = tree)

# sample from copula model
x = rHAC(1000, model)

# in the following case the true model is binary approximated
est.obj = estimate.copula(x, HAC_GUMBEL, ML, epsilon = 0)
plot(est.obj)

# consider also the aggregation of the variables
est.obj = estimate.copula(x, HAC_GUMBEL, ML, epsilon = 0.2)
plot(est.obj)
}
