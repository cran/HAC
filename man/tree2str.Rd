\name{tree2str}
\alias{tree2str}

\title{String structure of a HAC
}
\description{
The function prints the structure of HAC as string, so that the important characteristics of the copula can be identified.
}
\usage{
tree2str(hac, theta = TRUE, digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
  \item{theta}{
determines, whether the values of the dependency parameter(s) are printed (\code{TRUE}) or not (\code{FALSE}).
}
  \item{digits}{
  a positive integer value specifying the number of digits of the dependency parameter(s).
}
}
\value{
a string of the class \code{\link{character}}.
}
\author{
Ostap Okhrin
}

\seealso{
\code{\link{plot.hac}}
}
\examples{
# an object of the class hac is constructed
M = matrix(c(2, 0, 0, 0, 3, 4, 0, 0, ~X1 + X5, ~X2, ~X3, ~X4 + X6),
ncol = 4, byrow = TRUE)
g.model = hac(HAC_GUMBEL, M)

# the dependency parameters are returned within the curly brackets

# variables, which are hosted at the same node, are separated by a dot

# (X1.X5.X2)_{3} and (X3.X4.X6)_{4} are variables of the node with 
# dependency parameter 2
tree2str(g.model) # [1] "((X1.X5.X2)_{3}.(X3.X4.X6)_{4})_{2}"

# if theta = FALSE, only the structure of the variables is returned
tree2str(g.model, theta = FALSE) # [1] "((X1.X5.X2).(X3.X4.X6))"

# alternatively consider the following nested AC

# _{1} represents the initial node
# the first three variables are given by the subtree ((X4.X6.X2)_{3}.X3)_{2},
# X1 and X5
f.model = hac.full(type = HAC_CLAYTON, y = c(~X1 + X5, ~X3, ~X2, ~X4 + X6),
theta = c(1, 2, 3))
tree2str(f.model) # [1] "(((X4.X6.X2)_{3}.X3)_{2}.X1.X5)_{1}"
}