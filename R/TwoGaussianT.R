#' @name TwoGaussianT
#' @title A data generating model with three clusters 
#' @description  Generates data set consists of three clusters, two clusters from Gaussian distributions 
#' which are compact and closely located to each other but far from the third cluster generated from the Student's t 
#' distribution with wider spread.  
#' 
#' @usage TwoGaussianT(n)
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
#'
#'
#' @examples
#' dmat <- TwoGaussianT(1000)
#' plot2d(dmat$data, dmat$truelab)
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @importFrom stats rnorm
#' @export


TwoGaussianT <- function(n){
K <- 3
x1 <- rnorm(n/K, 0, .1)
y1 <- rnorm(n/K, 5, .1)

x2 <- rnorm(n/K, 0.5, .2)
y2 <- rnorm(n/K, 5.5, .2)

x3 <- rt(n/K, 25, 5)
y3 <- rt(n/K, 25, 10)

data <-cbind(c(x1,  x2, x3), c(y1, y2, y3))
truelab <- rep(c(1, 2, 3), c(n/K, n/K, n/K))
out <- list(n=nrow(data), K=K, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
