#' @name TChiGaussianF
#' 
#' @title A data generating model with five clusters  
#' @description  Generates data set consists of five clusters, one each from Gaussian, Student's t, Chi-squared, skew Gaussian,
#' and F distributions. 
#' @usage TChiGaussianF(n)
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
#'
#' @examples
#' dmat <- TChiGaussianF(500)
#' plot2d(dmat$data, dmat$truelab)
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' @importFrom sn rsn
#' @importFrom stats rnorm rchisq rf rt
#' @export


TChiGaussianF <- function(n){
K <- 5
x1 <- rchisq(n/K, df = 7, ncp = 50)
y1 <- rchisq(n/K, df = 10, ncp = 80)

x3 <- rf(n/K, df1 = 2, df2 = 6, ncp = 4)
y3 <- rf(n/K, df1 = 5, df2 = 5, ncp = 4)

x8 <- rnorm(n/K, 100, 0.9)
y8 <- rnorm(n/K, 0, 0.9)

x9 <- rt(n/K, 40, 100)
y9 <- rt(n/K, 35, 150)

x11 <-  rsn(n/K, xi = 20, omega = 0.9, alpha = 2, tau = 4)
y11 <-  rsn(n/K, xi = 200, omega = 0.8, alpha = 3, tau = 6)

data <- data.frame(x=c(x1, x3, x8, x9, x11), y=c(y1,  y3, y8, y9, y11))
truelab <- rep(1:K, each = n/K)

out <- list(n=n, K=K, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
