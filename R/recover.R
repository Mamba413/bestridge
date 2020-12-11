#' recover the coefficients from screening
#'
#' This function recovers the original coefficients in
#' the \code{beta.all} of the \code{bess} object when screening is utilized in \code{bess} function.
#'
#' When screening is used, the active set updates are restricted on a chosen subset and the \code{beta.all}
#'  component of \code{bess} object contains only coefficients on this subset. This function helps recover
#'  the whole coeffients.
#'
#' @param object A "\code{bess}" object.
#' @param sparse Logical or NULL, specifying whether the coefficients should be
#' presented as sparse matrix or not.
#' @author Canhong Wen, Aijun Zhang, Shijie Quan, Liyuan Hu, Kangkang Jiang, Yanhang Zhang, Jin Zhu and Xueqin Wang.
#' @seealso \code{\link{bess}}, \code{\link{plot.bess}},
#' \code{\link{bess}}.
#' @references Wen, C., Zhang, A., Quan, S. and Wang, X. (2020). BeSS: An R
#' Package for Best Subset Selection in Linear, Logistic and Cox Proportional
#' Hazards Models, \emph{Journal of Statistical Software}, Vol. 94(4).
#' doi:10.18637/jss.v094.i04.
#' @examples
#'
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
#' lm.bss <- bess(x, y, method = "sequential")
#'
#' recover(lm.bss)
#' recover(lm.bss, best.model = FALSE)
#'
#'@export

recover<- function(object, sparse=TRUE){
  if(object$method == "sequential"){
    beta.all <- lapply(object$beta.all, list.beta, object, sparse)
  } else{
    beta.all = matrix(0, length(object$beta), ncol = ncol(object$beta.all))
    beta.all[object$screening_A, ] = object$beta.all
    if(sparse){
      beta.all <- Matrix(beta.all)
    }
  }
  return(beta.all)
}

list.beta <- function(beta.mat, object, sparse){
  beta.all <- matrix(0, nrow=length(object$beta), ncol =  ncol(beta.mat))
  beta.all[object$screening_A, ] = beta.mat[[1]]
  if(sparse){
    beta.all <- Matrix(beta.all)
  }
  return(beta.all)
}
