#' @name TwoGaussian
#' @title A data generating model with two clusters  
#' @description  Generates data set consists of two clusters, one each from Gaussian distributions. 
#' @usage TwoGaussian(n)
#'
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
#' @examples
#' dmat <- TwoGaussian(1000)
#' plot2d(dmat$data, dmat$truelab)
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom stats rnorm
#' @export


TwoGaussian <- function(n){
K <- 2
x1 <- rnorm(n/K, 0, .7)
y1 <- rnorm(n/K, 5, .7)

x2 <- rnorm(n/K, +2, .1)
y2 <- rnorm(n/K, 5, .1)

data <- cbind(x=c(x1, x2), y=c(y1, y2))
truelab <- rep(c(1, 2), c(n/2, n/2))
out <- list(n=n, K=2, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
