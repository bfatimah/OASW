#' @name  C14D2
#' @title A data generating model with 14-clusters
#' @description Generates a data set consists of 14-clusters, one each from the specified Gaussian distributions. 
#' 
#' @usage C14D2(n)
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
#' dmat <- C14D2(350)
#' plot2d(dmat$data, dmat$truelab)
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom stats rnorm
#' @export


C14D2 <- function(n){
K <- 14
p <- 2

x1 <- rnorm(n/K, 0, .5)
y1 <- rnorm(n/K, 2, .5)
x8 <- rnorm(n/K, 0, .5)
y8 <- rnorm(n/K, -2, .5)
x6 <- rnorm(n/K, 4, .1)
y6 <- rnorm(n/K, 2, .1)
x9 <- rnorm(n/K, 3, .1)
y9 <- rnorm(n/K, -2, .1)
x2 <- rnorm(n/K, 2, .1)
y2 <- rnorm(n/K, 2, .1)
x7 <- rnorm(n/K, -4, .1)
y7 <- rnorm(n/K, 2, .1)
x10 <- rnorm(n/K, -3, .1)
y10 <- rnorm(n/K, -2, .1)
x3 <- rnorm(n/K, -2, .1)
y3 <- rnorm(n/K, 2, .1)
x4 <- rnorm(n/K, 3, .1)
y4 <- rnorm(n/K, 2, 0.7)
x5 <- rnorm(n/K, -3, .1)
y5 <- rnorm(n/K, 2, .7)
x11 <- rnorm(n/K, 2, .1)
y11 <- rnorm(n/K, -2, 0.7)
x12 <- rnorm(n/K, -2, .1)
y12 <- rnorm(n/K, -2, .7)
x13 <- rnorm(n/K, 4, .1)
y13 <- rnorm(n/K, -2, 0.7)
x14 <- rnorm(n/K, -4, .1)
y14 <- rnorm(n/K, -2, .7)
 
data <- cbind(c(x1, x2, x3, x4, x5, x6,x7,x8,x9,x10,x11,x12,x13, x14), 
              c(y1, y2, y3, y4, y5,y6,y7,y8,y9,y10,y11,y12,y13,y14))
truelab <- rep(1:K, each = n/K)

out <- list(n=n, K=K, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
