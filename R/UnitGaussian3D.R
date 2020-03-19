#' @name  UnitGaussian3D
#' @title A data generating model with nine clusters in three dimensions
#' @description Generates data set consists of nine clusters each drawn from a Gaussian distribution. See details 
#' for complete discription of the model's structure. 
#' 
#' @usage UnitGaussian3D(n)
#' @param n 
#' @return Returns a list having two components described as follows:
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has three iid dimensions.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#'
#' @examples
#' dmat <- UnitGaussian3D(891)
#' plot2d(dmat$data, dmat$truelab)
#' scatter3D(data[,1], data[,2], data[,3], colkey = FALSE, colvar = truelab, 
#'  ticktype = "detailed", pch = 16, bty = "b")
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @importFrom stats rnorm
#' @import plot3D
#' @export

UnitGaussian3D <- function(n){
    K <- 9
    p <- 3
    num <- (n/K)*3
    x <- matrix(rnorm(num, 0, 1),ncol=3)
    y <- x/sqrt(rowSums(x^2))
    
    x7 <- rnorm(n/K, -7, 0.1)
    y7 <- rnorm(n/K, -0.2, 0.1)
    z7 <- rnorm(n/K, -0.2, 0.1)
    C7 <- cbind(x7, y7, z7)
    
    x8 <- rnorm(n/K, -5.5, 0.6)
    y8 <- rnorm(n/K, 2.5, 0.8)
    z8 <- rnorm(n/K, 2.5, 0.6)
    C8 <- cbind(x8, y8, z8)
    
    x1 <- rnorm(n/K, -4, 0.4)
    y1 <- rnorm(n/K, -2.5, 0.3)
    z1 <- rnorm(n/K, -2.5, 0.4)
    C1 <- cbind(x1, y1, z1)
    
    x2 <- rnorm(n/K, 0.2, 0.1)
    y2 <- rnorm(n/K, -4, 0.1)
    z2 <- rnorm(n/K, -4, 0.1)
    C2 <- cbind(x2, y2, z2)
    
    x4 <- rnorm(n/K, 0.5, 0.1)
    y4 <- rnorm(n/K, 3, 0.1)
    z4 <- rnorm(n/K, 3, 0.1)
    C4 <- cbind(x4, y4, z4)
    
    x5 <- rnorm(n/K, 4.5, 0.6)
    y5 <- rnorm(n/K, -3, 0.8)
    z5 <- rnorm(n/K, -3, 0.6)
    C5 <- cbind(x5, y5, z5)
    
    x6 <- rnorm(n/K, 5, 0.4)
    y6 <- rnorm(n/K, 1.5, 0.3)
    z6 <- rnorm(n/K, 1.5, 0.4)
    C6 <- cbind(x6, y6, z6)
    
    x3 <- rnorm(n/K, 7, 0.1)
    y3 <- rnorm(n/K, -1, 0.1)
    z3 <- rnorm(n/K, -1, 0.1)
    C3 <- cbind(x3, y3, z3)
    
    data <- rbind(C1,  C2, C3, C4, C5, C6, C7, C8, y ) #y is not n/K its 99/3 or (n/K)/3
    #truelab <- rep(1:K, each = n/K)
    truelab <- rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9), c(n/K,n/K,n/K,n/K,n/K,n/K,n/K,n/K,n/K))
 
out <- list(n=nrow(data), K=K, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
