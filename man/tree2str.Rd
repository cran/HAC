\name{tree2str}
\alias{tree2str}

\title{String structure of HAC
}
\description{
The function prints the structure of HAC as string, so that the important characteristics of the copula can be identified.
}
\usage{
tree2str(hac, theta = TRUE, digits = 2)
}
\arguments{
  \item{hac}{
an object of the class \env{hac}.
}
  \item{theta}{
boolean. Determines, whether the values of the dependency parameter(s) are printed (\code{TRUE}) or not (\code{FALSE}).
}
  \item{digits}{
  a non-negative integer value specifying the number of digits of the dependency parameter(s).
}
}
\value{
a string of the class \code{\link{character}}.
}
\seealso{
\code{\link{plot.hac}}
}
\examples{
# a hac object is constructed
tree = list(list("X1", "X5", "X2", 3), list("X3", "X4", "X6", 4), 2)
model = hac(tree = tree)

# the dependency parameters are returned within the curly brackets
# variables nested at the same node are separated by a dot

tree2str(model) # [1] "((X1.X5.X2)_{3}.(X3.X4.X6)_{4})_{2}"

# (X1.X5.X2)_{3} and (X3.X4.X6)_{4} are the two variables nested at the
# initial node with dependency parameter 2

tree2str(model, theta = FALSE) # [1] "((X1.X5.X2).(X3.X4.X6))"

# if theta = FALSE, only the structure of the variables is returned

# alternatively consider the following nested AC

tree = list("X1", list("X5", "X2", 3), list("X3", "X4", "X6", 4), 1.01)
model = hac(tree = tree)

tree2str(model) # [1] "(X1.(X5.X2)_{3}.(X3.X4.X6)_{4})_{1.01}"

# _{1.01} represents the initial node
# the first three variables are given by the subtrees (X3.X4.X6)_{4},
# (X5.X2)_{3} and X1
}