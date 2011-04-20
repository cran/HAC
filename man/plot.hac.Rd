\name{plot.hac}
\alias{plot.hac}
\title{Plot of a HAC
}
\description{
The function plots the structure of (Hierarchical) Archimedean Copulae.
}
\usage{ 
\method{plot}{hac}(x, xlim, ylim, xlab, ylab, axes = FALSE,
col = "black", fg = "black", bg = "white", lwd = 2, index = FALSE, 
theta = TRUE, h = 0.45, l = 1.2, circles = 0.25, digits = 2, \dots)
}
\arguments{
  \item{x}{
a \env{hac} object. It can be constructed by \code{\link{hac}} or be the result of the estimation procedure (\code{\link{estimate.copula}}).
}
  \item{xlim, ylim}{
numeric vectors of the length 2, giving the limits of the x and y axes. The default values adjust the size of the coordinate plane automatically with respect to the dimension of the HAC. If the type of the \env{hac} is \code{AC_GUMBEL} or \code{AC_CLAYTON}, the limits of the y-axis have to be adjusted manually, e.g. \code{ylim = c(-0.4, 1.25)}.
}
  \item{xlab, ylab}{
titles for the x and y axes.
}
  \item{axes}{
the default value \code{FALSE} implies, that no axes appear.
}
  \item{col}{
defines the color of the lines, which connect the circles and rectangles.
}
  \item{fg}{
defines the color of the lines of the rectangles and circles equivalent to the color settings in \env{R}.
}
  \item{bg}{
defines the background color of the rectangles and circles equivalent to the color settings in \env{R}.
}
  \item{lwd}{
the width of the lines.
}
  \item{index}{
is boolean, i.e. \code{index = T} or \code{index = F}. If it is set to \code{TRUE}, a string, which represents the structure of the hosted variables, appears as subscript of the dependence parameter. If the dimension of the copula is bigger or equal than 8, the argument \code{index = TRUE} is ignored, i.e. \code{index = FALSE} by default.
}
  \item{theta}{
determines, whether the dependency parameter \eqn{\theta} or Kendall's rank correlation coefficient \eqn{\tau} is printed. 
}
  \item{h}{
the height of the rectangles.
}
  \item{l}{
the width of the rectangles.
}
  \item{circles}{
a positive number giving the radius of the circles.
}
  \item{digits}{
an integer specifying the number of digits of the dependence parameter.
}
  \item{\dots}{
arguments to be passed to methods, such as graphical parameters (see \code{\link{par}}).
}
}
\details{
The user is advised to make the plot window as big as possible, since some printed text (e.g. the correlation parameter) could be outside of the rectangles or circles. Subscripts can only appear, if the dimension of the copula is strictly smaller than 8. If the \env{hac} refers to Archimedean copulae, the variables are always labeled by ``V\eqn{i}'', where \eqn{i = 1, \dots, k} and \eqn{k} denotes the amount of variables. This is reasoned by the exchangeability of the variables in the case of AC.
}
\author{
Ostap Okhrin, Alexander Ristig
}
\seealso{
\code{\link{estimate.copula}}
}
\examples{
# a hac object is created
Matrix = matrix(c(1.988, 0, 0, 0, 3.234, 4.567, 0, 0, ~X1, ~X2, ~X3, ~X4),
nrow = 3, byrow = TRUE)
Object = hac(HAC_GUMBEL, Matrix)

# do not forget to make the plot window as big as possible
plot(Object, col = "red", fg = "red", index = TRUE)

# the same procedure works for an estimated object
x = rHAC(1000, Object)
est.obj = estimate.copula(x, HAC_GUMBEL, method = ML)

# the result can be properly presented without any effort
plot(est.obj, index = TRUE, digits = 5)

# as explained, ylim should be adjusted for AC
ac_model = hac(type = AC_GUMBEL, X = 3, dim = 3)
plot(ac_model, ylim = c(-0.4, 1.25), l = 1, h = 0.35)
}