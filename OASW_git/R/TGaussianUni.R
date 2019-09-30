#' @title TGaussianUni
#'
#' @description  Four clsuters two of which are generated from t, and Uniform, each and two from
#'  Gaussian distributions. The dimensions are i.i.d.
#'
#' @param  \code{n} number of observations in the dataset
#' @return Returns a list having two components described as follows:
#' \describe{
#'  \item{data}{dataset generated from a data generating process having \code{n} observations and \code{K} number of clusters}
#' \item{truelab}{clustering label vector corresponding to the known data generating process}
#' }
#' @details each cluster has equal number of observations. The dataset has 2 dimensions and K=4 number of clusters
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- TGaussianUni(100)
#' plot(dmat$data, col = dmat$truelab, pch = 16)
#' @export

TGaussianUni <- function(n){
    K <- 4
    x1 <- rt(n/K, 7, 10)
    y1 <- rt(n/K, 7, 30)
    
    x2 <- rnorm(n/K, 2, 1)
    y2 <- rnorm(n/K, 2, 1)
    
    x3 <- runif(n/K, 10, 15)
    y3 <- runif(n/K, 10, 15)
    
    x4 <- rnorm(n/K, 20, 0.1)
    y4 <- rnorm(n/K, 80, 2)
    
    data <- cbind(x=c(x1, x2, x3, x4), y=c(y1, y2, y3, y4))
    truelab <- rep(c(1, 2, 3, 4), each = n/K)

out <- list(data=data, truelab=truelab)
return(out)
}

