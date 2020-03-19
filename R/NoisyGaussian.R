#' @name NoisyGaussian
#' @title A data generating model with two clusters and added Uniform noise points
#' @description Generates data set consists of two clossely located Gaussian clusters with 500 Uniform noise points. 
#' 
#' @usage NoisyGaussian(n,noise.points)
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @param noise.points number of noie points requried to be added in data set
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has two iid dimensions. 
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @examples
#' dmat <- NoisyGaussian(200, 50)
#' plot(dmat$data, col=c("blue","green", "red")[dmat$truelab])
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom stats rnorm
#' @export


NoisyGaussian <- function(n, noise.points){
    K <- 2
    p <- 2
    
    x1 <- rnorm(n/K, 0, .1)
    y1 <- rnorm(n/K, 5, .1)
    
    x2 <- rnorm(n/K, 0.5, .2)
    y2 <- rnorm(n/K, 5.5, .2)
    
    x3 <- runif(noise.points, -10, 10)
    y3 <- runif(noise.points, -5, 15)
    
    data <- cbind(c(x1,  x2, x3), c(y1, y2, y3))
    truelab <- rep(c(1, 2, 3), c(n/K, n/K, noise.points))
    
out <- list(n=nrow(data), K=K, cluster_size=c(rep(n/K, each=K), noise.points), data=data, truelab=truelab)
return(out)
}
 
