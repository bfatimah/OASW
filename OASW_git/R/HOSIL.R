#' @name HOSIL
#' 
#' @title Hierarchical optimum average SILhouette width clustering
#'
#' @description HOSIL is a clustering alogromative hierarchical clustering algorithm.
#'  The cluster mearge are defined at each hierarchy using a new linkage method defined by
#'   the optimization of the ASW index. The method can also estimate the number of clusters 
#'   based on OASW linkage criterion. If number of clsuters is to be fixed
#'    (see \link{\code{fix.k}}) the minimum allowd is 3 and can be at most \code{(n-1)}.
#'
#' @usage HOSIL(dys,  distmethod = "euclidean", fix.K = 3)
#' 
#' @param dys A vector of pairwise distances between observations. Usually an object of class \code{"dist"} or a data matrix or data frame. In latter case also needs distance method. The default is set at \code{euclidean}. 
#' @param distmethod  the distance mehtod to be used. Current available methods are \code{euclidean}, \code{maximum}, \code{manhattan}, \code{canberra}, \code{binary}, or \code{minkowski}. See \code{\link[stats:dist]{dist}} for more details on these methods.
#' @param fix.K user defined number of clusters against which clustering is requried.
#' @param ... other parameters
#'
#' @return Returns a list with th efollwoing components
#' 
#'\describe{
#' \item{est_K}{estimated number of clusters }
#' \item{clus_vect_est}{clustering label vector for the estimated number of clusters}
#' \item{OASW_est_K}{Value of the objective function against the estimated number of clusters}
#' \item{all_OASW}{Values of the objective function (OASW linkage) against \code{n-1} to \code{2} number of clusters. Thus the first value \code{all_OASW[1]} gives the OASW linkage value against \code{(n-1)} number of clusters, the next value contains \code{(n-2)} then \code{(n-3),..., 3, 2}}
#' \item{all_clus}{A vector of length \code{n*(n-2)} containing the clustering label vectors for all number of clusters from \code{(n-1)} to 2, where n is number of observations. So first n values will be clustering labels for n-1 number of clusters the next n will be clustering lables for n-2 clusters and so on. For n number of observations if clustering for K number of clusters is requried do: \cr
#'  \code{r <- n-K} \cr
#'  \code{all_clus[(n*(r-1)+1):(n*r)] + 1 }\cr
#'  
#' to see clustering labels against that K. For the user convience fix.K is also created to reterive clustering against one of these K. } }
#'\describe{
#' If the number of clusters are specified by the user some additional output is also returned as follows
#' \item{fix_K}{user specified number of clusters }
#' \item{clus_vect_fix}{clustering label vector for the user specified number of clusters}
#' \item{OASW_fix_K}{Value of the objective function against the user specified number of clusters}
#'}
#'
#'
#' @author Fatima Batool \email{ ucakfba@ucl.ac.uk }
#'
#'
#' @examples
#' require(sn)
#' data_5 <- function(n){
#' K <- 6
#' x2 <- rexp(n/K, rate = 10)
#' y2 <- rexp(n/K, rate = 10)
#' shape1 <- 2
#' shape2 <- 3
#' x7 <- rbeta(n/K, shape1, shape2, ncp = 220)
#' y7 <- rbeta(n/K, shape1, shape2, ncp = 120)
#' shape = 10
#' x6 <- rweibull(n/K, shape, scale = 4)
#' y6 <- rweibull(n/K, shape, scale = 4)
#' rate  <-  2
#' x4 <- rgamma(n = n/K, shape = 15, rate )
#' y4 <- rgamma(n = n/K, shape = 15)
#' x10 <-runif(n/K, -6, -2)
#' y10 <- runif(n/K, -6, -2)
#' x11 <-  rsn(n/K, xi = 5, omega = 0.6, alpha = 4, tau = 5)
#' y11 <-  rsn(n/K,xi = 0, omega = 0.6, alpha = 4, tau = 5)
#' data <- data.frame(x=c(x2, x7, x6, x4, x11, x10), y=c(y2, y7, y6, y4, y11, y10))
#' truelab <- rep(c(1, 2, 3, 4, 5, 6), each = n/K)
#' out <- list(data=data, truelab=truelab)
#' return(out)
#' }
#' dmat <- data_5(96)$data
#' dys <- dist(dmat)
#' HOSIL(dys) 
#' oasw_clus <- HOSIL(dys, fix.K=6)
#' plot(dmat, col = oasw_clus$clus_vect_fix, pch = 16)
#' @references
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' Kaufman, L. and P. J. Rousseeuw (1990). Finding groups in data: an introduction to cluster analysis,
#'  Volume 344. John Wiley & Sons.
#'  
#' @export
HOSIL <- function(dys,  distmethod = "euclidean", fix.K = "NA", ...){
    
    if(inherits(dys, "dist") == FALSE){
        dys = dist(dys, method = distmethod)
    }
    
    n <- (1 + sqrt(1+8*length(dys)))/2
    K = n
    clus_lab = integer(n)
    best_clus_lab = integer(n)
    all_best_avg_silh = numeric(n-2)
    all_best_clus_lab =   integer(n*(n-2))
    clus_size = integer(K)
    dys_j = numeric(K)
    avg_dys_clus=  numeric(n*K)
    sil= numeric(n)
    copy_best_clus_lab =  numeric(n)
    disty <- filldys(dys)
    hosil_lab_swap(n, disty, clus_lab, clus_size,  dys_j, avg_dys_clus, sil, best_clus_lab, all_best_clus_lab, all_best_avg_silh, copy_best_clus_lab);
    
    est <- n - which.max(all_best_avg_silh)
    rk <- n - est
    clus_lab_est <-  all_best_clus_lab[(n*(rk-1)+1):(n*rk)] + 1
    oasw.value.est.k <- all_best_avg_silh[rk]
    est_k <- est

    #results for user defined K
   if(fix.K != "NA"){
            r <- n-fix.K
            clus_lab_fix <-  all_best_clus_lab[(n*(r-1)+1):(n*r)] + 1
            oasw.value.fix.K <- all_best_avg_silh[r]
        output_more <- list(fix.K=fix.K, OASW_fix_K=oasw.value.fix.K, clus_vect_fix=clus_lab_fix, est_K=est_k, clus_vect_est=clus_lab_est,  OASW_est_K=oasw.value.est.k,   all_OASW = all_best_avg_silh, all_clus=all_best_clus_lab)
       return(output_more)
    } else {
    output <- list(est_K=est_k,  OASW_est_K=oasw.value.est.k, clus_vect_est=clus_lab_est,  all_OASW = all_best_avg_silh, all_clus=all_best_clus_lab)
     return(output)
    }

}


 

