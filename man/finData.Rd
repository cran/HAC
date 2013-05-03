\name{finData}      
\docType{data}
\alias{finData}

\title{Financial data}

\description{
This data set consists the standardized residuals of the filtered daily log-returns of four oil corporations: Chevron Corporation (\code{CVX}), Exxon Mobil Corporation (\code{XOM}), Royal Dutch Shell (\code{RDSA}) and Total (\code{FP}), covering \eqn{n = 283} observations from 20110202 to 20120319. Intertemporal dependence is removed by usual ARMA-GARCH models, whose standardized residuals are used as \code{finData}.
}

\format{
A matrix containing 283 observations on 4 stocks. The tickers of the stocks are presented as \code{colnames}.
}

\source{Yahoo! Finance}

\examples{
# load the data
data(finData)
}
