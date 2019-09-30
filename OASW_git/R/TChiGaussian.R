#' @name TChiGaussian
#' 
#' @title Data model with five clusters  
#' @description  Must have a description. add some
#'
#'
#' @param  \code{n} number of observations in the dataset. These are equally divided between clusters
#' @return Returns a list having two components described as follows:
#' \describe{
#'  \item{data}{dataset generated from a data generating process having n observations and K number of clusters}
#' \item{truelab}{clustering label vector correspounding to the known data generating process}
#' }
#' @details each cluster has equal number of observations. The dataset has 2 dimensions and K=5 number of clusters
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- TChiGaussian(100)
#' plot(dmat$data, col = dmat$truelab, pch = 16)
#' @export

#require(sn)
TChiGaussian <- function(n){
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

out <- list(data=data, truelab=truelab)
return(out)
}
