#' @name ThreeMicroarrayMulti
#' @title A data generating model to simulate microarray-like settings 
#' @description  Generates data set consists of three clusters as defined in Van der Laan (2003).
#' 
#' @usage ThreeMicroarrayMulti(n)
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has 1000-dimensions and \code{K=7} number of clusters. The dimensions are iid.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#
#' @examples
#' dmat <- ThreeMicroarrayMulti(60)
#  plot2d(dmat$data, dmat$truelab)
#' @references 
#' Van der Laan, M., K. Pollard, and J. Bryan (2003). A new partitioning around medoids algorithm. 
#' \emph{Journal of Statistical Computation and Simulation}, 73(8), 575 584.
#' @import MASS
#' @export


ThreeMicroarrayMulti <- function(n){
    #see KP_dataset file for labels. this data also needs a transpose before clustering.
    #n <- 60 #n is for samples
    K <- 7 #K is number of clusters
    p <- 500 #p denote dimensions here 
    
    mean1 <- mean2 <- mean3 <- numeric(p)
    
    #First sub-population
    mean_1 <- rep(c(log(3), -log(3), 0), c(25, 25, 450))
    
    #Second sub-population
    mean_2 <- rep(c(0, log(3), -log(3), 0), c(50, 25, 25, 400))
    
    #Third sub-population
    mean_3 <- rep(c(0, log(3), -log(3), 0), c(100, 25, 25, 350))
    
    #covariance matrix
    sigma <- diag(log(1.6)^2, p, p)
    
    subpop1 <- mvrnorm(20, mean_1, sigma)
    subpop2 <- mvrnorm(20, mean_2, sigma)
    subpop3 <- mvrnorm(20, mean_3, sigma)
    
    data <- rbind(subpop1, subpop2, subpop3)
    truelab <- rep(c(1, 2, 3, 4, 5, 6, 7), c(25, 25, 25, 25, 25, 25, 350))
    
out <- list(n=n, K=2, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
