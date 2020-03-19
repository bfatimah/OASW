#' @name ElongUniGaussian
#' @title A data generating model with 9-clusters 
#' @description  Generates data set consists of 9-clusters, from Gaussian and Uniform distributions. See details for complete model definition. 
#' 
#' @usage ElongUniGaussian(n)
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has two dimensions.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- ElongUniGaussian(900)
#' plot2d(dmat$data, dmat$truelab)
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @importFrom stats rnorm
#' @export


ElongUniGaussian <- function(n){
    K <- 9
    p <- 2
    x1 <- rnorm(n/K, 0, 1)
    y1 <- rnorm(n/K, 8, 1)
    
    x2 <- runif(n/K, -2, 2)
    y2 <- runif(n/K, -2, 2)
    
    x3 <- rnorm(n/K, 0, 1)
    y3 <- rnorm(n/K, -8, 1)
    
    x4 <- rnorm(n/K, 8, 1)
    y4 <- rnorm(n/K, 0, 1)
    
    x5 <- rnorm(n/K, -8, 1)
    y5 <- rnorm(n/K, 0, 1)
    
    t <- runif(n/K, -1, 1) 
    x6 <- t + rnorm(n/K, -7, 0.2)
    y6 <- t + rnorm(n/K, -7, 0.2)
    
    x7 <- t +  rnorm(n/K, 7, 0.1) 
    y7 <- t +  rnorm(n/K, 7, 0.1) 
    
    t <- runif(n/K, -0.9, 0.9) 
    x8 <- t + rnorm(n/K, -7, 0.1)
    y8 <- t + rnorm(n/K, 7, 0.1)
    
    x9 <- t +  rnorm(n/K, 7, 0.1) 
    y9 <- t +  rnorm(n/K, -7, 0.1) 
    
    data <- cbind(c(x1, x2, x3, x4, x5, x6, x7, x8, x9), c(y1, y2, y3, y4, y5, y6, y7, y8, y9))
    truelab <- rep(1:K, each=n/K)
    
out <- list(n=nrow(data), K=K, cluster_size=c(rep(n/K, each=K)), data=data, truelab=truelab)
return(out)
}
 
