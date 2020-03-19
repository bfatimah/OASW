#' @name TNestedGaussian
#' @title  A data generating model with three nested clusters  
#' @description  Generates data set consists of three clusters. Two closely located Gaussain clusters are nested inside Student's t cluster.
#' 
#' @usage TNestedGaussian(n)
#' @param  n number of observations in the dataset. These are equally divided between clusters.
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
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @examples
#' dmat <- TNestedGaussian(300)
#' plot(dmat$data, col=c("blue", "red","green")[dmat$truelab], xlab = " ", ylab = " ")
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @importFrom stats rnorm
#' @export


TNestedGaussian <- function(n){
    K <- 3
    p <- 2
    
    x1 <- rnorm(n/K, 0, .1)
    y1 <- rnorm(n/K, 5, .1)
    
    x2 <- rnorm(n/K, 0.5, .2)
    y2 <- rnorm(n/K, 5.5, .2)
    
    x3 <- rt(n/K, 25, 0)
    y3 <- rt(n/K, 25, 5)
 
    data <- cbind(c(x1,  x2, x3), c(y1, y2, y3))
    truelab <- rep(c(1, 2, 3), c(n/K, n/K, n/K))
    
out <- list(n=nrow(data), K=K, cluster_size=c(n/K, n/K, n/K), data=data, truelab=truelab)
return(out)
}
 
