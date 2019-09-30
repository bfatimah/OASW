#' @name MultiDist
#' @title Six clusters model
#' @description  Must have a description. add some
#'
#'
#' @param  n number of observations in the dataset. These are equally divided between clusters
#' @return Returns a list having two components described as follows:
#' \describe{
#'  \item{data}{dataset generated from a data generating process having n observations and K number of clusters}
#' \item{truelab}{clustering label vector correspounding to the known data generating process}
#' }
#' @details each cluster has equal number of observations. The dataset has 2 dimensions and K=6 number of clusters
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- MultiDist(100)
#' plot(dmat$data, col = dmat$truelab, pch = 16)
#' @export

#library(sn)
MultiDist <- function(n){
    K <- 6
    x2 <- rexp(n/K, rate = 10)
    y2 <- rexp(n/K, rate = 10)
    shape1 <- 2
    shape2 <- 3
    x7 <- rbeta(n/K, shape1, shape2, ncp = 220)
    y7 <- rbeta(n/K, shape1, shape2, ncp = 120)
    shape = 10
    x6 <- rweibull(n/K, shape, scale = 4)
    y6 <- rweibull(n/K, shape, scale = 4)
    rate  <-  2
    x4 <- rgamma(n = n/K, shape = 15, rate )
    y4 <- rgamma(n = n/K, shape = 15)
    x10 <-runif(n/K, -6, -2)
    y10 <- runif(n/K, -6, -2)
    
    x11 <-  rsn(n/K, xi = 5, omega = 0.6, alpha = 4, tau = 5)
    y11 <-  rsn(n/K,xi = 0, omega = 0.6, alpha = 4, tau = 5)
    data <- data.frame(x=c(x2, x7, x6, x4, x11, x10), y=c(y2, y7, y6, y4, y11, y10))
    truelab <- rep(c(1, 2, 3, 4, 5, 6), each = n/K)
    
    out <- list(data=data, truelab=truelab)
    return(out)
}

