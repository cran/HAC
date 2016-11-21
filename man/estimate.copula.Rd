\name{estimate.copula}
\alias{estimate.copula}

\title{Estimation of Hierarchical Archimedean Copulae
}
\description{
The function estimates the parameters and determines the structure of Hierarchical Archimedean Copulae.
}
\usage{
estimate.copula(X, type = 1, method = 1, hac = NULL, epsilon = 0, 
agg.method = "mean", margins = NULL, na.rm = FALSE, max.min = TRUE, \dots)
}
\arguments{
  \item{X}{
a \eqn{n \times d} matrix. If there are no \code{colnames} provided, the names \code{X1, X2, ...} will be given.
}
  \item{type}{
  defines the copula family, see \code{\link{phi}} for an overview of implemented families.
}
  \item{method}{
the estimation method. Select between quasi Maximum Likelihood \code{1}, full Maximum Likelihood \code{2}, recursive Maximum Likelihood \code{3} and penalized Maximum Likelihood \code{4}.
}  
   \item{hac}{
a \code{\link{hac}} object, which determines the structure and provides initial values. An object must be provided, if \code{method = 2} referring to the full Maximum Likelihood procedure.
}
  \item{epsilon}{
  scalar \eqn{\geq 0}. The variables of consecutive nodes are aggregated, if the difference of the dependency parameters is smaller than \code{epsilon}. For a detailed explanation see also \code{\link{aggregate.hac}}.
}
  \item{agg.method}{
if \eqn{\epsilon > 0}, the new dependency parameter can be determined by \code{"mean"}, \code{"min"} or \code{"max"} of the two parameters, see \code{\link{aggregate.hac}}.
}
  \item{margins}{
  specifies the margins. The data matrix is assumed to contain the values of the marginal distributions by default, i.e. \code{margins = NULL}. If raw data are used, the margins can be determined nonparametrically, \code{"edf"}, or in a parametric way, e.g. \code{"norm"}. Following the latter approach, the parameters of the distributions are estimated by Maximum Likelihood. Building on these estimates the values of the univariate margins are computed. If the argument is defined as scalar, all margins are computed according to this specification. Otherwise, different margins can be defined, e.g. \code{c("norm", "t", "edf")} for a 3-dimensional sample. Almost all continuous functions of \code{\link{Distributions}} are available. Inappropriate usage of this argument might lead to misspecified margins, e.g. application of \code{"exp"} even though the sample contains negative values.
}
  \item{na.rm}{
boolean. If \code{na.rm = TRUE}, missing values, \code{NA}, contained in \code{X} are removed.
}
  \item{max.min}{
boolean. If \code{max.min = TRUE} and an element of \code{X} is \eqn{\geq 1} or \eqn{\leq 0}, it is set to \eqn{1-10^{-8}} and \eqn{10^{-8}} respectively.
}
  \item{\dots}{
further arguments passed to or from other methods, e.g. \code{\link{na.omit}}.
}
}
\value{
A \env{hac} object is returned.
}
\references{
Genest, C., Ghoudi, K., and Rivest, L. P. 1995, A Semiparametric Estimation Procedure of Dependence Parameters in Multivariate Families of Distributions, \emph{Biometrika} 82, 543-552.

Gorecki, J., Hofert, M. and Holena, M. 2014, On the Consistency of an Estimator for Hierarchical Archimedean Copulas, In Talaysova, J., Stoklasa, J., Talaysek, T. (Eds.) \emph{32nd International Conference on Mathematical Methods in Economics, Olomouc: Palacky University}, 239-244.

Joe, H. 2005, Asymptotic Efficiency of the Two-Stage Estimation Method for Copula-Based Models, \emph{Journal of Multivariate Analysis} 94(2), 401-419.

Okhrin, O., Okhrin, Y. and Schmid, W. 2013, On the Structure and Estimation of Hierarchical Archimedean Copulas, \emph{Journal of Econometrics} 173, 189-204.

Okhrin, O. and Ristig, A. 2014, Hierarchical Archimedean Copulae: The \verb{HAC} Package", \emph{Journal of Statistical Software}, 58(4), 1-20, \url{https://www.jstatsoft.org/v58/i04/}.

Okhrin, O., Ristig, A., Sheen J. and Trueck, S. 2015, Conditional Systemic Risk with Penalized Copula, \emph{SFB 649 Discussion Paper 2015-038, Sonderforschungsbereich 649, Humboldt University, Germany}.
}
\examples{
# define the copula model
tree = list(list("X1", "X5", 3), list("X2", "X3", "X4", 4), 2)
model = hac(type = 1, tree = tree)

# sample from copula model
x = rHAC(1000, model)

# in the following case the true model is binary approximated
est.obj = estimate.copula(x, type = 1, method = 1, epsilon = 0)
plot(est.obj)

# consider also the aggregation of the variables
est.obj = estimate.copula(x, type = 1, method = 1, epsilon = 0.2)
plot(est.obj)

# full ML estimation to yield more precise parameter
est.obj.full = estimate.copula(x, type = 1, method = 2, hac = est.obj)

# recursive ML estimation leads to almost identical results
est.obj.r = estimate.copula(x, type = 1, method = 3)
}

