#' @name ThreeMicroarray
#' @title  A data generating model with three microarray like clusters 
#' @description  Generates data set consists of three clusters. 
#' 
#' @usage ThreeMicroarray(n)
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has 1000-dimensions. The dimensions are iid.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @examples
#' dmat <- ThreeMicroarray(120)
#' plot2d(dmat$data, dmat$truelab)
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @import MASS
#' @export


ThreeMicroarray <- function(n){
    p <- 1000
    K <- 3
    mean1 <- mean2 <- mean3 <- numeric(p)
    #only mean of first 100 dimensions shifted to -3
    for(i in 1:100){
        mean1[i] <- -3
    }
    #may be needed if mean needs to shift from 0.
    for(i in 1:100){
        mean2[i] <- 0
    }
    
    for(i in 1:100){
        mean3[i] <- 3
    }
    
    sigma <- matrix(0, p, p)
    for(i in 1:p)
    {
        for(hh in 1:p){
            if(i == hh){
                sigma[i , hh] <- 1
            }
        }
    }
    
    subpop1 <- mvrnorm(n/K, mean1, sigma)
    subpop2 <- mvrnorm(n/K, mean2, sigma)
    subpop3 <- mvrnorm(n/K, mean3, sigma)
    
    data <- rbind(subpop1, subpop2, subpop3)
    truelab <- rep(1:K, each = n/K)
    
out <- list(n=n, K=K, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
