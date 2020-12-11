#' Best subset selection/Best subset ridge regression with a
#' specified model size and a shrinkage parameter
#'
#' Best subset selection with a specified model size for generalized
#' linear models and Cox's proportional hazard model.
#'
#'  Given a model size \eqn{s}, we consider the following best subset selection problem:
#'\deqn{\min_\beta -2 logL(\beta) ;{ s.t.} \|\beta\|_0 = s.}
#'And given a model size \eqn{s} and a shrinkage parameter \eqn{\lambda}
#', consider the following best subset ridge regression problem:
#'\deqn{\min_\beta -2 logL(\beta) + \lambda \Vert\beta \Vert_2^2; { s.t.} \|\beta\|_0 = s.}
#'
#'In the GLM case, \eqn{logL(\beta)} is the log-likelihood function;
#' In the Cox model, \eqn{logL(\beta)} is the log parital likelihood function.
#'
#'The best subset selection problem is solved by the primal dual active set algorithm,
#'see Wen et al. (2017) for details. This algorithm utilizes an active set updating strategy
#'via primal and dual vairables and fits the sub-model by exploiting the fact that their
#' support set are non-overlap and complementary.
#'
#' @param x Input matrix, of dimension \eqn{n \times p}; each row is an observation
#' vector and each column is a predictor/feature/variable.
#' @param y The response variable, of \code{n} observations. For \code{family = "binomial"} should be
#' a factor with two levels. For \code{family="poisson"}, \code{y} should be a vector with positive integer.
#'  For \code{family = "cox"}, \code{y} should be a two-column matrix
#' with columns named \code{time} and \code{status}.
#' @param type One of the two types of problems.
#' \code{type = "bss"} for the best subset selection,
#' and \code{type = "bsrr"} for the best subset ridge regression.
#' @param family One of the following models: \code{"gaussian"}, \code{"binomial"},
#' \code{"poisson"}, or \code{"cox"}. Depending on the response.
#' @param method The method to be used to select the optimal model size and \eqn{L_2} shrinkage. For
#' \code{method = "sequential"}, we solve the best subset selection and the best subset ridge regression
#' problem for each \code{s} in \code{1,2,...,s.max} and \eqn{\lambda} in \code{lambda.list}. For \code{method =
#' "gsection"}, which is only valid for \code{type = "bss"},
#' we solve the best subset selection problem with model size ranged between s.min and s.max,
#' where the specific model size to be considered is determined by golden section. we
#' solve the best subset selection problem with a range of non-continuous model
#' sizes. For \code{method = "pgsection"} and \code{"psequential"}, the Powell method is used to
#' solve the best subset ridge regression problem.
#' @param tune The criterion for choosing the model size and \eqn{L_2} shrinkage
#' parameters. Available options are \code{"gic"}, \code{"ebic"}, \code{"bic"}, \code{"aic"} and \code{"cv"}.
#' Default is \code{"gic"}.
#' @param s A specified model size
#' @param lambda A shrinkage parameter for \code{"bsrr"}.
#' @param screening.num Users can pre-exclude some irrelevant variables according to maximum marginal likelihood estimators before fitting a
#' model by passing an integer to \code{screening.num} and the sure independence screening will choose a set of variables of this size.
#' Then the active set updates are restricted on this subset.
#' @param normalize Options for normalization. \code{normalize = 0} for
#' no normalization. Setting \code{normalize = 1} will
#' only subtract the mean of columns of \code{x}.
#' \code{normalize = 2} for scaling the columns of \code{x} to have \eqn{\sqrt n} norm.
#' \code{normalize = 3} for subtracting the means of the columns of \code{x} and \code{y}, and also
#' normalizing the columns of \code{x} to have \eqn{\sqrt n} norm.
#' If \code{normalize = NULL}, by default, \code{normalize} will be set \code{1} for \code{"gaussian"},
#' \code{2} for \code{"binomial"}, \code{3} for \code{"cox"}.
#' @param weight Observation weights. Default is \code{1} for each observation.
#' @param max.iter The maximum number of iterations in the bess function.
#' In most of the case, only a few steps can guarantee the convergence. Default
#' is \code{20}.
#' @param warm.start Whether to use the last solution as a warm start. Default
#' is \code{TRUE}.
#' @param nfolds The number of folds in cross-validation. Default is \code{5}.
#' @param group.index A vector of integers indicating the which group each variable is in.
#' For variables in the same group, they should be located in adjacent columns of \code{x}
#' and their coorespinding index in \code{group.index} should be the same.
#' Denote the first group as \code{1}, the second \code{2}, etc.
#' If you do not fit a model with a group structure,
#' please set \code{group.index = NULL}. Default is \code{NULL}.
#' @param seed seed to be used to devide the sample into K cross-validation folds. Default is \code{NULL}.
#' @return A list with class attribute 'bess' and named components:
#' \item{beta}{The best fitting coefficients.} \item{coef0}{The best fitting
#' intercept.}
#' \item{loss}{The training loss of the fitting model.}
#' \item{s}{The model size.}
#' \item{lambda}{The shrinkage parameter.}
#' \item{family}{Type of the model.}
#' \item{nsample}{The sample size.}
#' \item{type}{Either \code{"bss"} or \code{"bsrr"}.}
#' @author Canhong Wen, Aijun Zhang, Shijie Quan, Liyuan Hu, Kangkang Jiang, Yanhang Zhang, Jin Zhu and Xueqin Wang.
#' @seealso \code{\link{plot.bess}}, \code{\link{summary.bess}}, \code{\link{recover.beta}},
#' \code{\link{coef.bess}}, \code{\link{predict.bess}}.
#' @references Wen, C., Zhang, A., Quan, S. and Wang, X. (2020). BeSS: An R
#' Package for Best Subset Selection in Linear, Logistic and Cox Proportional
#' Hazards Models, \emph{Journal of Statistical Software}, Vol. 94(4).
#' doi:10.18637/jss.v094.i04.
#' @examples
#'
#' #-------------------linear model----------------------#
#' # Generate simulated data
#' n <- 200
#' p <- 20
#' k <- 5
#' rho <- 0.4
#' SNR <- 10
#' cortype <- 1
#' seed <- 10
#' Data <- gen.data(n, p, k, rho, family = "gaussian", cortype = cortype, snr = SNR, seed = seed)
#' x <- Data$x[1:140, ]
#' y <- Data$y[1:140]
#' x_new <- Data$x[141:200, ]
#' y_new <- Data$y[141:200]
#' lm.bss <- bess.one(x, y, s = 5)
#' lm.bsrr <- bess.one(x, y, type = "bsrr", s = 5, lambda = 0.01)
#' coef(lm.bss)
#' coef(lm.bsrr)
#' print(lm.bss)
#' print(lm.bsrr)
#' summary(lm.bss)
#' summary(lm.bsrr)
#' pred.bss <- predict(lm.bss, newx = x_new)
#' pred.bsrr <- predict(lm.bsrr, newx = x_new)
#'
#' #-------------------logistic model----------------------#
#' #Generate simulated data
#' Data = gen.data(n, p, k, rho, family = "binomial", cortype = cortype, snr = SNR, seed = seed)
#'
#' x <- Data$x[1:140, ]
#' y <- Data$y[1:140]
#' x_new <- Data$x[141:200, ]
#' y_new <- Data$y[141:200]
#' logi.bss <- bess.one(x, y, family = "binomial", s = 5)
#' logi.bsrr <- bess(x, y, type = "bsrr", family = "binomial", s = 5, lambda = 0.01)
#' coef(logi.bss)
#' coef(logi.bsrr)
#' print(logi.bss)
#' print(logi.bsrr)
#' summary(logi.bss)
#' summary(logi.bsrr)
#' pred.bss <- predict(logi.bss, newx = x_new)
#' pred.bsrr <- predict(logi.bsrr, newx = x_new)
#'
#' #-------------------coxph model----------------------#
#' #Generate simulated data
#' Data <- gen.data(n, p, k, rho, family = "cox", scal = 10)
#'
#' x <- Data$x[1:140, ]
#' y <- Data$y[1:140, ]
#' x_new <- Data$x[141:200, ]
#' y_new <- Data$y[141:200, ]
#' cox.bss <- bess(x, y, family = "cox", s = 5)
#' cox.bsrr <- bess(x, y, type = "bsrr", family = "cox", s = 5, lambda = 0.01)
#' coef(cox.bss)
#' coef(cox.bsrr)
#' print(cox.bss)
#' print(cox.bsrr)
#' summary(cox.bss)
#' summary(cox.bsrr)
#' pred.bss <- predict(cox.bss, newx = x_new)
#' pred.bsrr <- predict(cox.bsrr, newx = x_new)
#'#----------------------High dimensional linear models--------------------#
#` p <- 1000
#` data <- gen.data(n, p, family = "gaussian", cortype = cortype, snr = SNR, seed = seed)
#`
#` # Best subset selection with SIS screening
#` fit <- bess.one(data$x, data$y, s= 10, screening.num = 100)
#`
#'#-------------------group selection----------------------#
#'beta <- rep(c(rep(1,2),rep(0,3)), 4)
#'Data <- gen.data(n, p, rho=0.4, beta = beta, snr = 100, seed =10)
#'x <- Data$x
#'y <- Data$y
#'
#'group.index <- c(rep(1, 2), rep(2, 3), rep(3, 2), rep(4, 3),
#'                 rep(5, 2), rep(6, 3), rep(7, 2), rep(8, 3))
#'lm.group <- bess(x, y, s.min=1, s.max = 8, type = "bss", group.index = group.index, s = 5)
#'lm.groupbsrr <- bess(x, y, type = "bsrr", s = 5, lambda = 0.01, group.index = group.index)
#'coef(lm.group)
#'coef(lm.groupbsrr)
#'print(lm.group)
#'print(lm.groupbsrr)
#'#'summary(lm.group)
#'summary(lm.groupbsrr)
#'pred.group <- predict(lm.group, newx = x_new)
#'pred.groupl0l2 <- predict(lm.groupbsrr, newx = x_new)
#'
#' @export


bess.one <- function(x, y, family = c("gaussian", "binomial", "poisson", "cox"), type = c("bss", "bsrr"),
                             s, lambda= 0,
                             screening.num = NULL,
                             normalize = NULL, weight = NULL,
                             max.iter = 20, warm.start = TRUE,
                             nfolds = 5,
                             group.index =NULL,
                             seed=NULL){
  family <- match.arg(family)
  type <- match.arg(type)

  res <- bess(x, y, family = family, type = type,
                       method ="sequential",
                       tune = "gic",
                       s.list=s, lambda.list = lambda,
                       s.min=s, s.max=s,
                       lambda.min = lambda, lambda.max = lambda, nlambda = 1,
                       screening.num = screening.num,
                       normalize = normalize, weight = weight,
                       max.iter = max.iter, warm.start = warm.start,
                       nfolds = nfolds,
                       group.index =group.index,
                       seed=seed)
  res$s <- s
  res$bess.one <- TRUE
  res$call <- match.call()
  return(res)
}
