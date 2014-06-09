\name{to.logLik}
\alias{to.logLik}

\title{log-likelihood
}
\description{
\code{to.logLik} returns either the log-likehood function depending on a vector \code{theta} for a given sample \code{X} or the value of the log-likelihood, if \code{eval = TRUE}.
}
\usage{
to.logLik(X, hac, eval = FALSE, margins = NULL, sum.log = TRUE, 
na.rm = FALSE, \dots)
}
\arguments{
  \item{X}{
a data matrix. The number of columns and the corresponding names have to coincide with the specifications of the copula model \code{hac}. The sample \code{X} has to contain at least \eqn{2} rows (observations), as the values of the underlying density cannot be computed otherwise.
} 
  \item{hac}{
an object of the class \code{\link{hac}}.
}
  \item{eval}{
boolean. If \code{eval = FALSE}, the non-evaluated log-likelihood function depending on a parameter vector \code{theta} is returned and one default argument, the \code{density}, is returned. The values of \code{theta} are increasingly ordered.
}
 \item{margins}{
  specifies the margins. The data matrix \code{X} is assumed to contain the values of the marginal distributions by default, i.e. \code{margins = NULL}. If raw data are used, the margins can be determined nonparametrically, \code{"edf"}, or in parametric way, e.g. \code{"norm"}. See \code{\link{estimate.copula}} for a detailed explanation.
}
 \item{sum.log}{
boolean. If \code{sum.log = FALSE}, the values of the individual log-likelihood contributions are returned.
}
  \item{na.rm}{
boolean. If \code{na.rm = TRUE}, missing values, \code{NA}, contained in \code{X} are removed.
}
  \item{\dots}{
arguments to be passed to \code{\link{na.omit}}.
}
}
\seealso{
\code{\link{dHAC}}
}
\examples{
# construct a hac-model
tree = list(list("X1", "X5", 3), list("X2", "X3", "X4", 4), 2)
model = hac(type = 1, tree = tree)

# sample from copula model
sample = rHAC(1000, model)

# check the accurancy of the estimation procedure
ll = to.logLik(sample, model)
ll.value = to.logLik(sample, model, eval = TRUE)

ll(c(2, 3, 4)) == ll.value # [1] TRUE
}
