% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deviance.bsrr.R
\name{deviance.bsrr}
\alias{deviance.bsrr}
\title{Extract the deviance from a "bsrr.one" object.}
\usage{
\method{deviance}{bsrr}(object, best.model = TRUE, ...)
}
\arguments{
\item{object}{A "\code{bsrr}" object.}

\item{best.model}{Whether only return the loglikelihood of the best model. Default is \code{TRUE}.
If \code{best.model = FALSE}, the loglikelihood of the best models with model size and
 \eqn{\lambda} in the original \code{s.list} and \code{lambda.list} (for \code{method = "sequential"})
 or in the iteration path (for \code{method = "gsection"}, \code{method = "pgsection"},
 and \code{method = "psequential"}) is returned.}

\item{\dots}{additional arguments}
}
\value{
A matrix or vector containing the deviance for each model is returned.
For \code{bsrr} object fitted by \code{sequantial} method, values in each row in the
returned matrix corresponding to the model size in \code{s.list}, and each column the shrinkage parameters
in \code{lambda.list}.

For \code{bsrr} object fitted by \code{gsection}, \code{pgsection} and \code{psequential}, the returned vector
contains deviance for fitted models in each iteration. The coefficients of those model can be extracted
from \code{beta.all} and \code{coef0.all} in the \code{bsrr} object.
}
\description{
Similar to other deviance methods, which returns deviance from a fitted "\code{bsrr.one}" object.
}
\examples{

# Generate simulated data
n <- 200
p <- 20
k <- 5
rho <- 0.4
seed <- 10
Tbeta <- rep(0, p)
Tbeta[1:k*floor(p/k):floor(p/k)] <- rep(1, k)
Data <- gen.data(n, p, k, rho, family = "gaussian", seed = seed)
lm.bsrr <- bsrr(Data$x, Data$y, method = "sequential")

deviance(lm.bsrr)
deviance(lm.bsrr, best.model = FALSE)
}
\seealso{
\code{\link{bsrr}}, \code{\link{summary.bsrr}}.
}
\author{
Liyuan Hu, Kangkang Jiang, Yanhang Zhang, Jin Zhu, Canhong Wen and Xueqin Wang.
}
