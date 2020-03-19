#' @name TenNest
#' @title A data generating model with ten nested clusters  
#' @description Generates data set consists of ten clusters, one each from Gaussian distributions. See details for full model definition.
#' 
#' @usage TenNest(n)
#' @param n number of observations in the dataset. These are equally divided between clusters
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has five thousands iid dimensions. 
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @examples
#' dmat <- TenNest(1000)
#' plot2d(dmat$data, dmat$truelab)
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @importFrom stats rnorm
#' @export


TenNest <- function(n){
    p <- 5000
    K <- 10
    n1 <- n/K
    data <- matrix(numeric(n*p),  ncol = p)
    meansvec <- sample(c(-21, -18, -15, -9, -6, 6, 9, 15, 18, 21), 10)
    meansvec <- sample(c(-16, -13, -10, -6, -3, 3, 6, 10, 13, 16), 10)
    mu1 <- meansvec[1] 
    mu2 <- meansvec[2] 
    mu3 <- meansvec[3] 
    mu4 <- meansvec[4] 
    mu5 <- meansvec[5] 
    mu6 <- meansvec[6] 
    mu7 <- meansvec[7] 
    mu8 <- meansvec[8] 
    mu9 <- meansvec[9] 
    mu10 <- meansvec[10] 
    
    x1 <- rnorm(n1, mu1, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))	
    x2 <- rnorm(n1, mu2, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x3 <- rnorm(n1, mu3, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x4 <- rnorm(n1, mu4, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x5 <- rnorm(n1, mu5, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x6 <- rnorm(n1, mu6, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x7 <- rnorm(n1, mu7, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x8 <- rnorm(n1, mu8, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x9 <- rnorm(n1, mu9, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    x10 <- rnorm(n1, mu10, sample(c(0.005, 0.1, 0.2,  0.3, 0.4), 1))
    
    datainit <- c(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
    #datainit <- c(x3, x9, x7, x10, x6, x2, x8, x1, x4, x5)
    
    for(v in 1:p){
        for(i in 1:n){
            data[i, v] <- datainit[i]
        }
    }
    
truelab <- rep(1:10, each = n1)
out <- list(n=n, K=2, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
