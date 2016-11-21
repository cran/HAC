\name{dHAC, pHAC, rHAC}
\alias{dHAC}
\alias{pHAC}
\alias{rHAC}
\title{pdf, cdf and random sampling
}
\description{
\code{dHAC} and \code{pHAC} compute the values of the copula's density and cumulative distribution function respectively. \code{rHAC} samples from HAC.
}
\usage{
dHAC(X, hac, eval = TRUE, margins = NULL, na.rm = FALSE, \dots)
pHAC(X, hac, margins = NULL, na.rm = FALSE, \dots)
rHAC(n, hac)
}
\arguments{
  \item{X}{
a data matrix. The number of columns and the corresponding names have to coincide with the specifications of the copula model \code{hac}.
} 
  \item{hac}{
an object of the class \code{\link{hac}}.
}
  \item{n}{
number of observations.
}
  \item{margins}{
  specifies the margins. The data matrix \code{X} is assumed to contain the values of the marginal distributions by default, i.e. \code{margins = NULL}. If raw data are used, the margins can be determined nonparametrically, \code{"edf"}, or in parametric way, e.g. \code{"norm"}. See \code{\link{estimate.copula}} for a detailed explanation.
}

  \item{na.rm}{
boolean. If \code{na.rm = TRUE}, missing values, \code{NA}, contained in \code{X} are removed.
}
  \item{eval}{
boolean. If \code{eval = FALSE}, a non-evaluated \code{\link{function}} is returned. Note, that \code{\link{attr}} \code{"gradient"} of the returned function corresponds to the values density.
}
  \item{\dots}{
arguments to be passed to \code{\link{na.omit}}.
}
}

\value{
    \code{rHAC} retruns a \eqn{n \times d} matrix, where \eqn{d} refers to the dimension of the HAC. \code{dHAC} and \code{pHAC} return vectors. The computation of the density might be time consuming for high-dimensions, since the density is defined as \eqn{d}-th derivative of the HAC with respect to its arguments \eqn{u_1, \ldots, u_d}. 
}
\details{
    Sampling schemes of hierarchical and densities of simple Archimedean copula are based on functions of the \env{copula} package.
}
\references{
Hofert, M. 2011, Efficiently Sampling Nested Archimedean Copulas, \emph{Computational Statistics & Data Analysis} 55, 57-70.

Joe, H. 1997, Multivariate Models and Dependence Concepts, \emph{Chapman & Hall}.

McNeil, A. J. 2008, Sampling Nested Archimedean Copulas, \emph{Journal of Statistical Computation and Simulation} 78, 567-581.

Nelsen, R. B. 2006, An Introduction to Copulas, \emph{Spinger}, 2nd Edition.

Okhrin, O. and Ristig, A. 2014, Hierarchical Archimedean Copulae: The \verb{HAC} Package", \emph{Journal of Statistical Software}, 58(4), 1-20, \url{https://www.jstatsoft.org/v58/i04/}.

Savu, C. and Trede, M. 2010, Hierarchies of Archimedean copulas, \emph{Quantitative Finance} 10, 295-304.
}
\seealso{
\code{\link{estimate.copula}}, \code{\link{to.logLik}}
}
\examples{
# AC example
# define the underlying model
model = hac(type = 4, tree = list("X1", "X2", 2))

# sample from model
sample = rHAC(100, model)

# returns the pdf/cdf at each vector of the sample
d.values = dHAC(sample, model)
p.values = pHAC(sample, model)

# HAC example
# the underlying model
y = c("X1", "X2", "X3")
theta = c(1.5, 3)
model = hac.full(type = 1, y, theta)

# define sample from copula model
sample = rHAC(100, model)

# returns the pdf/cdf at each point of the sample
d.values = dHAC(sample, model)
p.values = pHAC(sample, model)

# construct a hac-model
tree = list(list("X1", "X5", 3), list("X2", "X3", "X4", 4), 2)
model = hac(type = 1, tree = tree)

# sample from copula model
sample = rHAC(1000, model)

# check the accurancy of the estimation procedure
result1 = estimate.copula(sample)
result2 = estimate.copula(sample, epsilon = 0.2)
}

