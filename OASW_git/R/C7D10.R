#' @name C7D10
#' @title Seven clusters in ten dimensions having unequal within cluster variations.
#' @description  All the clusters are generated from Gaussain distributions. The dimensions are i.i.d.
#'
#'
#' @param  \code{n} number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows:
#' \describe{
#'  \item{data}{dataset generated from a data generating process having \code{n} observations and \code{K} number of clusters}.
#' \item{truelab}{clustering label vector correspounding to the known data generating process}
#' }
#' @details each cluster has equal number of observations. The dataset has 10 dimensions and \code{K=7} number of clusters.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- C7D10(175)
#' plot(dmat$data, col = dmat$truelab, pch = 16)
#' @export


C7D10 <- function(n){
K <- 7
x1 <- rnorm(n/K, 0, .3)
y1 <- rnorm(n/K, 5, .3)
C1 <- cbind(x1, y1, y1+3, y1+6, y1+9, y1+12, y1+15, y1+18, y1+21, y1+24)

x4 <- rnorm(n/K, 0, .2)
y4 <- rnorm(n/K, 3.5, .5)
C4 <- cbind(x4, y4, y4+3, y4+6, y4+9, y4+12, y4+15, y4+18, y4+21, y4+24)

x5 <- rnorm(n/K, -0.5, .1)
y5 <- rnorm(n/K,3.5, .3)
C5 <- cbind(x5, y5, y5+3, y5+6, y5+9, y5+12, y5+15, y5+18, y5+21, y5+24)

x6 <- rnorm(n/K, 0.5, .2)
y6 <- rnorm(n/K,3.5, .3)
C6 <- cbind(x6, y6, y6-3, y6-6, y6-9, y6-12, y6-15, y6-18, y6-21, y6-24)

x7 <- rnorm(n/K, 0, .2)
y7 <- rnorm(n/K,6.5, .5)
C7 <- cbind(x7, y7, y7+3, y7+3, y7+3, y7+3, y7-3, y7-3, y7-3, y7-3)

x8 <- rnorm(n/K, -0.5, .1)
y8 <- rnorm(n/K,6.5, .3)
C8 <- cbind(x8, y8, y8-3, y8-6, y8-9, y8-12, y8-15, y8-18, y8-21, y8-24)

x9 <- rnorm(n/K, 0.5, .2)
y9 <- rnorm(n/K,6.5, .3)
C9 <- cbind(x9, y9, y9-3, y9-6, y9-9, y9-12, y9-15, y9-18, y9-21, y9-24)
truelab <- rep(1:K, each = n/K)
data <- rbind(C1, C4, C5, C6, C7, C8, C9)

out <- list(data=data, truelab=truelab)
return(out)

}
