\name{aggregate.hac}
\alias{aggregate.hac}

\title{Aggregation of variables
}
\description{
\code{aggregate} tests, whether the absolute difference of the parameters of two subsequent nodes is smaller than a constant, i.e. \eqn{\vert \theta_{2} - \theta_{1} \vert < \epsilon}, where \eqn{\theta_{i}} denotes the dependency parameter with \eqn{\theta_{2} < \theta_{1}}, \eqn{\epsilon \geq 0}. If the absolute difference is smaller than the constant, the variables of the nodes are aggregated in a single node with new dependency parameter, e.g. \eqn{ \theta_{new} = (\theta_{1} + \theta_{2})/2}. This procedure is applied to all consecutive nodes of the HAC \code{x}.
}
\usage{
\method{aggregate}{hac}(x, epsilon = 0, method = "mean", \dots)
}
\arguments{
  \item{x}{
an object of the class \env{hac}.
}
  \item{epsilon}{
scalar \eqn{\geq 0}.
}
  \item{method}{
determines, whether the new parameter is the \code{"mean"}, \code{"min"} or \code{"max"} of the fused parameters.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\value{
an object of the class \env{hac}.
}
\seealso{
\code{\link{hac}}
}
\examples{
# Example 1:
# an object of the class hac is constructed, whose parameters are close
copula = hac(type = 1, tree = list("X1", list("X2", "X3", 2.05), 2))

# the function aggregate returns a simple Archimedean copula

copula_ag = aggregate(copula, epsilon = 0.1)
tree2str(copula_ag) # [1] "(X1.X2.X3)_{2.02}"

# the structure does not change for a smaller epsilon

copula_ag = aggregate(copula, epsilon = 0.01)
tree2str(copula_ag) # [1] "((X2.X3)_{2.05}.X1)_{2}"

# Example 2:
# consider the binary tree

Object = hac.full(type = 1, y = c("X1", "X2", "X3", "X4", "X5"), 
theta = c(1.01, 1.02, 2, 2.01))

tree2str(Object) # [1] "((((X5.X4)_{2.01}.X3)_{2}.X2)_{1.02}.X1)_{1.01}"

# applying aggregate.hac with epsilon = 0.02 leads to

Object_ag = aggregate(Object, 0.02)
tree2str(Object_ag) # [1] "((X3.X5.X4)_{2}.X1.X2)_{1.02}"
}