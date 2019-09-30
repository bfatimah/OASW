#'  @name ThreeGaussian
#'  @title ThreeGaussian
#'
#' @description  Three Gaussian clusters of unequal within cluster variations
#'
#' @param  \code{n} number of observations in the dataset. These are equally divided between clusters
#' @return Returns a list having two components described as follows:
#' \describe{
#'  \item{data}{dataset generated from a data generating process having \code{n} observations and \code{K} number of clusters}
#' \item{truelab}{clustering label vector correspounding to the known data generating process}
#' }
#' @details each cluster has equal number of observations. The dataset has 2 dimensions and K=3 number of clusters
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- ThreeGaussian(100)
#' plot(dmat$data, col = dmat$truelab, pch = 16)
#' @export


ThreeGaussian <- function(n){
    K <- 3
    x1 <- rnorm(n/K, 0, .7)
    y1 <- rnorm(n/K, 0, .7)
    
    x2 <- rnorm(n/K, -2, .1)
    y2 <- rnorm(n/K, 0, .1)
    
    x3 <- rnorm(n/K, 2, .1)
    y3 <- rnorm(n/K, 0, .1)
    
    data <- cbind(x=c(x1, x2, x3), y=c(y1, y2, y3))
    truelab <- rep(c(1, 2, 3), c(n/K, n/K, n/K))
    
    out <- list(data = data, truelab=truelab)
    return(out)
}

