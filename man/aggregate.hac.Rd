\name{aggregate.hac}
\alias{aggregate.hac}

\title{Aggregation of variables
}
\description{
The function \code{Aggregate} tests, whether the absolute difference of the dependency parameters of two successive nodes is smaller than a constant, i.e. \eqn{\vert \theta_{1} - \theta_{2} \vert < \epsilon}, where \eqn{\theta_{i}} denote the dependency parameter with \eqn{\theta_{1} < \theta_{2}}, \eqn{\epsilon \geq 0} and \eqn{i \in \{1, 2 \}}. If the difference is smaller than the constant, the variables of the nodes are aggregated in one node of the HAC with a new dependency parameter \eqn{(\theta_{1} + \theta_{2}) / 2}. \code{Aggregate} applies this procedure to all consecutive nodes of a HAC.
}
\usage{
\method{aggregate}{hac}(x, epsilon = 0.01, \dots)
}
\arguments{
  \item{x}{
an object of the class \env{hac}.
}
  \item{epsilon}{
a scalar \eqn{\geq 0}.
}
  \item{\dots}{
arguments to be passed to methods.
}
}
\details{
The flexibility of HAC can lead to relatively complicated structures, if the amount of variables is high. In these cases it is possible, that the function does not aggregate variables, although it is theoretically possible to aggregate them. This is reasoned by the binary approximation of the HAC, i.e. a node can host at most 2 successive nodes. Such a situation is illustrated in example 3.
}
\value{
an object of the class \env{hac}.
}

\author{
Alexander Ristig
}
\seealso{
\code{\link{hac}}
}
\examples{
# Example 1:
# an object of the class hac is constructed, such that
# the values of the dependency parameters are closely

M = matrix(c(1.01, 0, 0, 0, ~X1 + X2, 1.02, 0, 0, 0, 0, ~X3 + X4, ~X5 + X6), 
ncol = 4, byrow = TRUE)
Object = hac(HAC_GUMBEL, M)

# the function aggregate returns an AC_GUMBEL model, since the tree
# collapses to an Archimedean copulae

Object_ag = aggregate(Object, 0.1)

# if we choose epsilon smaller, e.g. 0.0001, nothing occurs and the
# model remains as it is

Object_not_ag = aggregate(Object, 0.0001)

# Example 2:
# consider a binary tree

Object = hac.full(y = c(~X1, ~X2, ~X3, ~X4, ~X5), theta = c(1.01, 1.02, 2, 2.01))
tree2str(Object) # [1] "((((X5.X4)_{2.01}.X3)_{2}.X2)_{1.02}.X1)_{1.01}"

# after the application of aggregate.hac, with epsilon = 0.011, two nodes remain

Object_ag = aggregate(Object, 0.011)
tree2str(Object_ag) # [1] "((X5.X4.X3)_{2.005}.X2.X1)_{1.015}"

# Example 3:
# consider the following deeply nested model
X = matrix(c(0.7, 0, 0, 0, 0, 0, 0, 0,
1, 0.75, 0, 0, 0, 0, 0, 0,
1.01, 1.02, ~X1, ~X2, 0, 0, 0, 0,
~X3, 2, 2.01, 2.02, 0, 0, 0, 0,
0, 0, ~X4, ~X5, ~X6, ~X7, ~X8, ~X9), nrow = 5, byrow = TRUE)
Object = hac(HAC_CLAYTON, X)

# to get an idea about the structure, plot Object

plot(Object)

# if aggregate.hac is applied with epsilon = 0.1, the object simplifies
# but the nodes with the parameter values 1.01 and 1.02 remain, since
# the upper node cannot host 3 nodes and 1 variable

Object_ag = aggregate(Object, 0.1)
plot(Object_ag)
}