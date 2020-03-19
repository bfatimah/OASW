#' @name  CCo4D2
#' @title A data generating model with four clusters
#' @description Generates a data set consists of 4-clusters, one each from the specified Gaussian distributions. 
#' 
#' @usage CCo4D2(n)
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has two iid dimensions and \code{K=4} clusters. 
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- CCo4D2(1000)
#' plot2d(dmat$data, dmat$truelab)
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @importFrom stats rnorm
#' @export


CCo4D2 <- function(n){
K <- 4
x1 <- rnorm(n/K, 1, 0.2)
y1 <- rnorm(n/K, 1, 0.2)
x2 <- rnorm(n/K, 1.5, 0.2)
y2 <- rnorm(n/K, 1.5,  0.2)
x3 <- rnorm(n/K, 1.5, 0.2)
y3 <- rnorm(n/K, 1, 0.2)
x4 <- rnorm(n/K, 1, 0.2)
y4 <- rnorm(n/K, 1.5,  0.2)
 
data <- cbind(c(x1, x2, x3, x4), c(y1, y2, y3, y4))
truelab <- rep(1:K, each=n/K)
out <- list(n=n, K=2, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
