\name{phi, phi.inv}
\alias{phi}
\alias{phi.inv}

\title{Generator function
}
\description{The Archimedean generator function and its inverse.
}
\usage{
phi(x, theta, type)
phi.inv(x, theta, type)
}
\arguments{
  \item{x}{
a scalar, vector or matrix at which the function is evaluated. The support of the functions has to be taken into account, i.e. \eqn{x \in [0, \infty]} for the generator function and \eqn{x \in [0, 1]} for its inverse.
}
  \item{theta}{
  the feasible copula parameter, i.e. \eqn{\theta \in [1, \infty)} for the Gumbel and Joe family, \eqn{\theta \in (0, \infty)} for the Clayton and Frank family and \eqn{\theta \in [0, 1)} for the Ali-Mikhail-Haq family.
 
}
\item{type}{
select between the following integer numbers for specifying the \code{type} of the hierarchical Archimedean copula (HAC) or Archimedean copula (AC):
	\itemize{
        	\item \code{ 1} = HAC Gumbel
	        \item \code{ 2} = AC Gumbel
	        \item \code{ 3} = HAC Clayton
	        \item \code{ 4} = AC Clayton
        	\item \code{ 5} = HAC Frank
	        \item \code{ 6} = AC Frank
        	\item \code{ 7} = HAC Joe
	        \item \code{ 8} = AC Joe
        	\item \code{ 9} = HAC Ali-Mikhail-Haq
	        \item \code{10} = AC Ali-Mikhail-Haq
}}
}
\examples{
x = runif(100, min = 0, max = 100)
phi(x, theta = 1.2, type = 1)

# do not run
# phi(x, theta = 0.8, type = 1) 
# In phi(x, theta = 0.8, type = 1) : theta >= 1 is required.
}
