#' @name Co5d5
#' @title A data generating model with two clusters   
#' @description  Generates data set consists of five correlated clusters, one each from Gaussian distributions. 
#' 
#' @usage Co5d5(n)
#' @param  n number of observations in the dataset. These are equally divided between clusters.
#' @return Returns a list having two components described as follows;
#'  \describe{
#'  \item{n}{number of observations.}
#'   \item{K}{number of clusters.}
#'    \item{cluster size}{number of observations in each cluster.}
#'  \item{data}{data set generated from the model.}
#' \item{truelab}{clustering label vector correspounding to the known data generating model.}
#' }
#' @details The data set has five correlated dimensions. 
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- Co5d5(1000)
#' plot2d(dmat$data, dmat$truelab)
#' pairplots(dmat$data, dmat$truelab) 
#' 
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' @import MASS
#' @import ggplot2
#' @import GGally
#' @export

Co5d5 <- function(n){
    p <- 5
    K <- 5
    sgsq1 <- 9
    sgsq2 <- 17
    sgsq3 <- 12
    sgsq4 <- 2
    sgsq5 <- 16
    r12 <-  r21<- 0.1
    r13 <- r31 <- 0.1
    r14 <- r41 <- 0.1
    r15 <- r51 <- 0.1
    r23 <- r32 <- 0.1
    r24 <-  r42 <- 0.1
    r25 <-  r52 <- 0.1
    r34 <- r43 <- 0.1
    r35 <- r53 <- 0.1
    r45 <-  r54 <- 0.1
    sigma1 = matrix(c(sgsq1, r12*sqrt(sgsq1)*sqrt(sgsq2), r13*sqrt(sgsq1)*sqrt(sgsq3), r14*sqrt(sgsq1)*sqrt(sgsq4), r15*sqrt(sgsq1)*sqrt(sgsq5), 
                      r21*sqrt(sgsq1)*sqrt(sgsq2), sgsq2, r23*sqrt(sgsq2)*sqrt(sgsq3), r24*sqrt(sgsq2)*sqrt(sgsq4), r25*sqrt(sgsq2)*sqrt(sgsq5),
                      r31*sqrt(sgsq1)*sqrt(sgsq3), r32*sqrt(sgsq2)*sqrt(sgsq3), sgsq3, r34*sqrt(sgsq3)*sqrt(sgsq4), r35*sqrt(sgsq3)*sqrt(sgsq5),
                      r41*sqrt(sgsq1)*sqrt(sgsq4), r42*sqrt(sgsq2)*sqrt(sgsq4),  r34*sqrt(sgsq3)*sqrt(sgsq4), sgsq4, r45*sqrt(sgsq4)*sqrt(sgsq5),
                      r51*sqrt(sgsq1)*sqrt(sgsq5), r52*sqrt(sgsq2)*sqrt(sgsq5),  r53*sqrt(sgsq3)*sqrt(sgsq5), r54*sqrt(sgsq4)*sqrt(sgsq5), sgsq5 ), 5, 5)
    
    mean1 <- c(0, 0, 0, 0, 0)
    C1 <- mvrnorm(n/K, mean1, sigma1)
    
    sgsq1 <- 1
    sgsq2 <- 1
    sgsq3 <- 1
    sgsq4 <- 1
    sgsq5 <- 1
    r12 <-  r21<- 0.3
    r13 <- r31 <- 0.3
    r14 <- r41 <- -0.3
    r15 <- r51 <-  -0.3
    r23 <- r32 <--0.3
    r24 <-  r42 <- 0.3
    r25 <-  r52 <- -0.3
    r34 <- r43 <- 0.3
    r35 <- r53 <- -0.3
    r45 <-  r54 <- -0.3
    
    sigma2 = matrix(c(sgsq1, r12*sqrt(sgsq1)*sqrt(sgsq2), r13*sqrt(sgsq1)*sqrt(sgsq3), r14*sqrt(sgsq1)*sqrt(sgsq4), r15*sqrt(sgsq1)*sqrt(sgsq5), 
                      r21*sqrt(sgsq1)*sqrt(sgsq2), sgsq2, r23*sqrt(sgsq2)*sqrt(sgsq3), r24*sqrt(sgsq2)*sqrt(sgsq4), r25*sqrt(sgsq2)*sqrt(sgsq5),
                      r31*sqrt(sgsq1)*sqrt(sgsq3), r32*sqrt(sgsq2)*sqrt(sgsq3), sgsq3, r34*sqrt(sgsq3)*sqrt(sgsq4), r35*sqrt(sgsq3)*sqrt(sgsq5),
                      r41*sqrt(sgsq1)*sqrt(sgsq4), r42*sqrt(sgsq2)*sqrt(sgsq4),  r34*sqrt(sgsq3)*sqrt(sgsq4), sgsq4, r45*sqrt(sgsq4)*sqrt(sgsq5),
                      r51*sqrt(sgsq1)*sqrt(sgsq5), r52*sqrt(sgsq2)*sqrt(sgsq5),  r53*sqrt(sgsq3)*sqrt(sgsq5), r54*sqrt(sgsq4)*sqrt(sgsq5), sgsq5 ), 5, 5)
    
    mean2 <- c(5, 10, 3, 7, 6)
    C2 <- mvrnorm(n/K, mean2, sigma2)
    
    sgsq1 <- 25
    sgsq2 <- 9
    sgsq3 <- 16
    sgsq4 <- 1
    sgsq5 <- 49
    r12 <-  r21<- 0.2
    r13 <- r31 <- 0.2
    r14 <- r41 <- -0.2
    r15 <- r51 <-  -0.2
    r23 <- r32 <--0.2
    r24 <-  r42 <- 0.2
    r25 <-  r52 <- -0.2
    r34 <- r43 <- 0.2
    r35 <- r53 <- -0.2
    r45 <-  r54 <- -0.2
    sigma3 = matrix(c(sgsq1, r12*sqrt(sgsq1)*sqrt(sgsq2), r13*sqrt(sgsq1)*sqrt(sgsq3), r14*sqrt(sgsq1)*sqrt(sgsq4), r15*sqrt(sgsq1)*sqrt(sgsq5), 
                      r21*sqrt(sgsq1)*sqrt(sgsq2), sgsq2, r23*sqrt(sgsq2)*sqrt(sgsq3), r24*sqrt(sgsq2)*sqrt(sgsq4), r25*sqrt(sgsq2)*sqrt(sgsq5),
                      r31*sqrt(sgsq1)*sqrt(sgsq3), r32*sqrt(sgsq2)*sqrt(sgsq3), sgsq3, r34*sqrt(sgsq3)*sqrt(sgsq4), r35*sqrt(sgsq3)*sqrt(sgsq5),
                      r41*sqrt(sgsq1)*sqrt(sgsq4), r42*sqrt(sgsq2)*sqrt(sgsq4),  r34*sqrt(sgsq3)*sqrt(sgsq4), sgsq4, r45*sqrt(sgsq4)*sqrt(sgsq5),
                      r51*sqrt(sgsq1)*sqrt(sgsq5), r52*sqrt(sgsq2)*sqrt(sgsq5),  r53*sqrt(sgsq3)*sqrt(sgsq5), r54*sqrt(sgsq4)*sqrt(sgsq5), sgsq5 ), 5, 5)
    
    mean3 <- c(15, 70, 50, 55, 80)
    C3 <- mvrnorm(n/K, mean3, sigma3)
    
    sgsq1 <- 5
    sgsq2 <- 0.9
    sgsq3 <- 1.6
    sgsq4 <- 1
    sgsq5 <- 4.9
    r12 <-  r21<- 0.1
    r13 <- r31 <- 0.1
    r14 <- r41 <- -0.7
    r15 <- r51 <-  -0.2
    r23 <- r32 <--0.2
    r24 <-  r42 <- 0.2
    r25 <-  r52 <- -0.9
    r34 <- r43 <- 0.2
    r35 <- r53 <- -0.2
    r45 <-  r54 <- -0.2
    sigma4 = matrix(c(sgsq1, r12*sqrt(sgsq1)*sqrt(sgsq2), r13*sqrt(sgsq1)*sqrt(sgsq3), r14*sqrt(sgsq1)*sqrt(sgsq4), r15*sqrt(sgsq1)*sqrt(sgsq5), 
                      r21*sqrt(sgsq1)*sqrt(sgsq2), sgsq2, r23*sqrt(sgsq2)*sqrt(sgsq3), r24*sqrt(sgsq2)*sqrt(sgsq4), r25*sqrt(sgsq2)*sqrt(sgsq5),
                      r31*sqrt(sgsq1)*sqrt(sgsq3), r32*sqrt(sgsq2)*sqrt(sgsq3), sgsq3, r34*sqrt(sgsq3)*sqrt(sgsq4), r35*sqrt(sgsq3)*sqrt(sgsq5),
                      r41*sqrt(sgsq1)*sqrt(sgsq4), r42*sqrt(sgsq2)*sqrt(sgsq4),  r34*sqrt(sgsq3)*sqrt(sgsq4), sgsq4, r45*sqrt(sgsq4)*sqrt(sgsq5),
                      r51*sqrt(sgsq1)*sqrt(sgsq5), r52*sqrt(sgsq2)*sqrt(sgsq5),  r53*sqrt(sgsq3)*sqrt(sgsq5), r54*sqrt(sgsq4)*sqrt(sgsq5), sgsq5 ), 5, 5)
    
    mean4 <- c(70, 80, 70, 70, 70)
    C4 <- mvrnorm(n/K, mean4, sigma4)
    
    sgsq1 <- 2
    sgsq2 <- 9
    sgsq3 <- 3
    sgsq4 <- 1
    sgsq5 <- 4
    r12 <-  r21<- 0.2
    r13 <- r31 <- 0.2
    r14 <- r41 <- -0.3
    r15 <- r51 <-  -0.1
    r23 <- r32 <--0.1
    r24 <-  r42 <- 0.2
    r25 <-  r52 <- -0.1
    r34 <- r43 <- 0.1
    r35 <- r53 <- -0.2
    r45 <-  r54 <- -0.9
    
    sigma5 = matrix(c(sgsq1, r12*sqrt(sgsq1)*sqrt(sgsq2), r13*sqrt(sgsq1)*sqrt(sgsq3), r14*sqrt(sgsq1)*sqrt(sgsq4), r15*sqrt(sgsq1)*sqrt(sgsq5), 
                      r21*sqrt(sgsq1)*sqrt(sgsq2), sgsq2, r23*sqrt(sgsq2)*sqrt(sgsq3), r24*sqrt(sgsq2)*sqrt(sgsq4), r25*sqrt(sgsq2)*sqrt(sgsq5),
                      r31*sqrt(sgsq1)*sqrt(sgsq3), r32*sqrt(sgsq2)*sqrt(sgsq3), sgsq3, r34*sqrt(sgsq3)*sqrt(sgsq4), r35*sqrt(sgsq3)*sqrt(sgsq5),
                      r41*sqrt(sgsq1)*sqrt(sgsq4), r42*sqrt(sgsq2)*sqrt(sgsq4),  r34*sqrt(sgsq3)*sqrt(sgsq4), sgsq4, r45*sqrt(sgsq4)*sqrt(sgsq5),
                      r51*sqrt(sgsq1)*sqrt(sgsq5), r52*sqrt(sgsq2)*sqrt(sgsq5),  r53*sqrt(sgsq3)*sqrt(sgsq5), r54*sqrt(sgsq4)*sqrt(sgsq5), sgsq5 ), 5, 5)
    
    mean5 <- c(55, 55, 55, 55, 55)
    C5 <- mvrnorm(n/K, mean5, sigma5)
    data <- rbind(C1, C2, C3, C4, C5)
    truelab <- rep(1:K, each = n/K)
    
out <- list(n=n, K=K, cluster_size=n/K, data=data, truelab=truelab)
return(out)
}
 
