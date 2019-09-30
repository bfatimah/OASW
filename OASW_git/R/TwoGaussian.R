#' @name TwoGaussian
#' @title TwoGaussian
#' @description  Two Gaussian clusters of unequall within cluster variations.
#' 
#'@usage TwoGaussian(n=100)
#'
#' @param  \code{n} number of observations in the dataset. These are equally divided between clusters
#' @return Returns a list having two components described as follows:
#' \describe{
    #'  \item{data}{dataset generated from a data generating process having \code{n} observations and \code{K=2} number of clusters}
    #' \item{truelab}{clustering label vector correspounding to the known data generating process}
    #' }
#' @details each cluster has equal number of observations. The dataset has 2 dimensions and K=2 number of clusters
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- TwoGaussian(100)
#' plot(dmat$data, col = dmat$truelab, pch = 16)
#' @export


TwoGaussian <- function(n){
    
x1 <- rnorm(n/2, 0, .7)
y1 <- rnorm(n/2, 5, .7)

x2 <- rnorm(n/2, +2, .1)
y2 <- rnorm(n/2, 5, .1)

data <- cbind(x=c(x1, x2), y=c(y1, y2))
truelab <- rep(c(1, 2), c(n/2, n/2))
out <- list(data=data, truelab=truelab)
return(out)
}
 
