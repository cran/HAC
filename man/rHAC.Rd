\name{rHAC}
\alias{rHAC}

\title{Simulation of a defined HAC
}
\description{
The simulation of HAC is implemented in this function.
}
\usage{
rHAC(n, hac)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the number of generated random vectors.
}
  \item{hac}{
an object of the class \code{\link{hac}}.
}
}
\value{
A \eqn{n \times d} matrix, where \eqn{d} refers to the dimension of the HAC.
}
\references{
Hofert, M. 2011, Efficiently Sampling Nested Archimedean Copulas, \emph{Computational Statistics & Data Analysis} 55, 57-70.

McNeil, A.,J. 2008, Sampling Nested Archimedean Copulas, \emph{Journal of Statistical Computation and Simulation} 78, 567-581.
}
\author{
Ostap Okhrin, Alexander Ristig
}
\seealso{
\code{\link{estimate.copula}}
}
\examples{
ac_model = hac(type = AC_GUMBEL, X = 3, dim = 3)
sample = rHAC(1000, ac_model)
plot(ac_model, ylim = c(-0.4, 1.25), l = 1, h = 0.35)

M = matrix(c(2, 0, 0, 0, 3, 4, 0, 0, ~X1 + X5, ~X2, ~X3, ~X4 + X6),
ncol = 4, byrow = TRUE)
g.model = hac(HAC_CLAYTON, M)
plot(g.model, index = TRUE)

# sample from copula g.model
x = rHAC(5000, g.model)
}