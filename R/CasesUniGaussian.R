#' @name UniGauss
#' @title A data generating model with two clusters
#' @description Generates a data set consists of 2-clusters, one each from the Gaussian and Uniform distributions. 
#' @usage UniGauss(n)
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
#' @examples
#' dmat <- UniGauss(1000)
#' plot2d(dmat$data, dmat$truelab)
#' plot(dmat$data, col=dmat$truelab, xlim=c(-20, 20), ylim=c(-20, 20))
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom stats rnorm runif
#' @export

UniGauss <- function(n){
K <- 2
x1 <- rnorm(n/K, 0, 1)
y1 <- rnorm(n/K, 5, 1)

x2 <- runif(n/K, -10, 10)
y2 <- runif(n/K, -1, 1)

data <- cbind(c(x1, x2), c(y1, y2))
truelab <- rep(c(1, 2), c(n/K, n/K))

out <- list(n=nrow(data), K=K, cluster_size=c(n/K, n/K), data=data, truelab=truelab)
return(out)
}

#' @name ShortUniGauss
#' @title A data generating model with two clusters
#' @description Generates a data set consists of 2-clusters, one each from the Gaussian and Uniform distributions.   
#' @usage ShortUniGauss(n)
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
#' @examples
#' dmat <- ShortUniGauss(1000)
#' plot2d(dmat$data, dmat$truelab)
#' plot(dmat$data, col=dmat$truelab, xlim=c(-20, 20), ylim=c(-20, 20))
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom stats rnorm runif
#' @export

ShortUniGauss  <- function(n){
K <- 2
x1 <- rnorm(n/K, 0, 1)
y1 <- rnorm(n/K, 5, 1)
x2 <- runif(n/K, -5, 5)
y2 <- runif(n/K, -1, 1)
data <- cbind(c(x1, x2), c(y1, y2))
truelab <- rep(c(1, 2), c(n/K, n/K))

out <- list(n=nrow(data), K=K, cluster_size=c(n/K, n/K), data=data, truelab=truelab)
return(out)
}


#' @name LongUniGauss
#' @title A data generating model with two clusters
#' @description Generates a data set consists of 2-clusters, one each from the Gaussian and Uniform distributions. 
#' @usage LongUniGauss(n)
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
#' @examples
#' dmat <- LongUniGauss(1000)
#' plot2d(dmat$data, dmat$truelab)
#' plot(dmat$data, col=dmat$truelab, xlim=c(-20, 20), ylim=c(-20, 20))
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom stats rnorm runif
#' @export
#' 
LongUniGauss <- function(n){
K <- 2
x1 <- rnorm(n/K, 0, 1)
y1 <- rnorm(n/K, 5, 1)
x2 <- runif(n/K, -15, 15)
y2 <- runif(n/K, -1, 1)
data <- cbind(c(x1, x2), c(y1, y2))
truelab <- rep(c(1, 2), c(n/K, n/K))

out <- list(n=nrow(data), K=K, cluster_size=c(n/K, n/K), data=data, truelab=truelab)
return(out)
}

#' @name FarUniGauss
#' @title A data generating model with two clusters
#' @description Generates a data set consists of 2-clusters, one each from the Gaussian and Uniform distributions.  
#' @usage FarUniGauss(n)
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
#' @examples
#' dmat <- FarUniGauss(1000)
#' plot2d(dmat$data, dmat$truelab)
#' plot(dmat$data, col=dmat$truelab, xlim=c(-20, 20), ylim=c(-20, 20))
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom stats rnorm runif
#' @export

FarUniGauss <- function(n){
    K <- 2
    x1 <- rnorm(n/K, 0, 1)
    y1 <- rnorm(n/K, 10, 1)
    
    x2 <- runif(n/K, -10, 10)
    y2 <- runif(n/K, -1, 1)
    
    data <- cbind(c(x1, x2), c(y1, y2))
    truelab <- rep(c(1, 2), c(n/K, n/K))
    
    out <- list(n=nrow(data), K=K, cluster_size=c(n/K, n/K), data=data, truelab=truelab)
    return(out)
}
