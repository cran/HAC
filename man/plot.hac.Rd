\name{plot.hac}
\alias{plot.hac}
\title{Plot of a HAC
}
\description{
The function plots the structure of Hierarchical Archimedean Copulae.
}
\usage{ 
\method{plot}{hac}(x, xlim = NULL, ylim = NULL, xlab = "", ylab = "", 
col = "black", fg = "black", bg = "white", col.t = "black", lwd = 2, 
index = FALSE, numbering = FALSE, theta = TRUE, h = 0.4, l = 1.2, 
circles = 0.25, digits = 2, \dots)
}
\arguments{
  \item{x}{
a \env{hac} object. It can be constructed by \code{\link{hac}} or be the result of \code{\link{estimate.copula}}.
}
  \item{xlim, ylim}{
numeric vectors of length 2, giving the limits of the x and y axes. The default values adjust the size of the coordinate plane automatically with respect to the dimension of the HAC.
}
  \item{xlab, ylab}{
titles for the x and y axes.
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
  \item{col.t}{
defines the text color equivalent to the color settings in \env{R}.
}
  \item{lwd}{
the width of the lines.
}
  \item{index}{
boolean. If \code{index = TRUE}, strings, which illustrate the subcopulae of the nodes, are used as subsrcipts of the dependency parameters.
}
  \item{numbering}{
boolean. If \code{index = TRUE} and \code{numbering = TRUE}, the dependency parameters are numbered. If \code{x} is returned by \code{estimate.copula}, the numbers correpsond to the estimation stages. 
}
  \item{theta}{
boolean. Determines, whether the dependency parameter \eqn{\theta} or Kendall's rank correlation coefficient \eqn{\tau} is printed. 
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
arguments to be passed to methods, e.g. graphical parameters (see \code{\link{par}}).
}
}
\references{
Okhrin, O. and Ristig, A. 2014, Hierarchical Archimedean Copulae: The \verb{HAC} Package", \emph{Journal of Statistical Software}, 58(4), 1-20, \url{https://www.jstatsoft.org/v58/i04/}.
}
\seealso{
\code{\link{estimate.copula}}
}
\examples{
# a hac object is created

tree = list(list("X1", "X5", 3), list("X2", "X3", "X4", 4), 2)
model = hac(type = 1, tree = tree)
plot(model)

# the same procedure works for an estimated object

sample = rHAC(2000, model)
est.obj = estimate.copula(sample, epsilon = 0.2)
plot(est.obj)
}
