\name{rAC}
\alias{rAC}
\title{Simulation of Archimedean Copulae
}
\description{
The Marshal-Olkin Algorithm is implemented to simulate from Archimedean Copulae.
}
\usage{
rAC(n, theta = 1.5, dim = 2, type = AC_GUMBEL)
}

\arguments{
  \item{n}{
the number of generated random vectors.
}
  \item{theta}{
the dependency parameter of the AC.
}
  \item{dim}{
the dimension of the copula.
}
  \item{type}{
choose the AC-type between \code{AC_GUMBEL} and \code{AC_CLAYTON}.
}
}
\value{
A \eqn{n \times d} matrix, where d refers to the dimension of the AC.
}
\references{Marshall, A.W. and Olkin, I. 1988,
	Families of Multivariate Distributions,
	\emph{Journal of the American Statistical Association} 83, 834-841.
}
\author{
Ostap Okhrin
}
\seealso{
\code{\link{rHAC}}
}
\examples{
# a sample of the 3 dimensional Clayton copula
sample = rAC(5000, 2.5, 3, AC_CLAYTON)
}