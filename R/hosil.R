#' @name hosil
#' 
#' @title Hierarchical optimum average silhouette width clustering
#' @description hosil is a clustering alogromative hierarchical clustering algorithm.
#'  The cluster mearge are defined at each hierarchy using a new linkage method defined by
#'   the optimization of the ASW index. The method can also estimate the number of clusters 
#'   based on OASW linkage criterion. If number of clsuters is to be fixed
#'    (see \code{fixK}) the minimum allowed is 3 and can be at most \code{(n-1)}.
#'
#' @usage hosil(dys,  distmethod = "euclidean", fixK = "NA")
#' 
#' @param dys A vector of pairwise distances between observations. Usually an object of class 
#' \code{"dist"} or a data matrix or data frame. In latter case also needs distance method. 
#' The default is set at \code{euclidean}. 
#' @param distmethod  the distance mehtod to be used. Current available methods are
#'  \code{euclidean}, \code{maximum}, \code{manhattan}, \code{canberra}, \code{binary}, or 
#'  \code{minkowski}. See \code{\link[stats:dist]{dist}} for more details on these methods.
#' @param fixK user defined number of clusters against which clustering is requried.
# @param ... other parameters
#'
#' @return Returns a list with the follwoing components
#'\describe{
#' \item{est_K}{estimated number of clusters.}
#' \item{OASW_est_K}{Value of the objective function against the estimated number of clusters.}
#' \item{clus_vect_est}{clustering label vector for the estimated number of clusters.}
#' \item{all_OASW}{Values of the objective function (OASW linkage) against \code{n-1} to \code{2} 
#' number of clusters. Thus the first value \code{all_OASW[1]} gives the OASW linkage value 
#' against \code{(n-1)} number of clusters, the next value contains \code{(n-2)} 
#' then \code{(n-3),..., 3, 2}.}
# \item{all_clus}{A vector of length \code{n*(n-2)} containing the clustering label vectors
#  for all number of clusters from \code{(n-1)} to 2, 
#  where n is number of observations. First n values will be clustering labels for n-1 number 
#  of clusters the next n will be clustering lables for n-2 clusters and so on. 
#  For n number of observations if clustering for K number of clusters is requried do: \cr
#  \code{r <- n-K} \cr
#  \code{all_clus[(n*(r-1)+1):(n*r)] + 1} \cr
# to see clustering labels against that K. 
# For the user convience fix.K is also created to reterive clustering against one of these K. In addition,  If the number of clusters are specified by the user some additional output is also returned as follows: }
#'
#' \item{fix_K}{user specified number of clusters }
#' \item{clus_vect_fix}{clustering label vector for the user specified number of clusters}
#' \item{OASW_fix_K}{Value of the objective function against the user specified number of clusters}
#'}
#'
#' @author Fatima Batool \email{ ucakfba@ucl.ac.uk }
#'
#'
#' @examples
#' dmat <- MultiDist(350)
#' dys <- dist(dmat$data)
#' oasw_clustering <- hosil(dys) 
#' oasw_clustering <- hosil(dys, fixK=6)
#' plot(dmat, col = oasw_clus$clus_vect_fix, pch = 16)
#' @references
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' Kaufman, L. and P. J. Rousseeuw (1990). Finding groups in data: an introduction to cluster analysis,
#'  Volume 344. John Wiley & Sons.
#' 
#' @importFrom Rcpp sourceCpp evalCpp
#' @export

#-----------------
#hosil algorithm
#-----------------
hosil <- function(dys,  distmethod = "euclidean", fixK = "NA"){
    
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
   if(fixK != "NA"){
            r <- n-fixK
            clus_lab_fix <-  all_best_clus_lab[(n*(r-1)+1):(n*r)] + 1
            oasw.value.fix.K <- all_best_avg_silh[r]
        output_more <- list(fix_K=fixK, OASW_fix_K=oasw.value.fix.K, clus_vect_fix=clus_lab_fix, est_K=est_k, clus_vect_est=clus_lab_est,  OASW_est_K=oasw.value.est.k,   all_OASW = all_best_avg_silh, all_clus=all_best_clus_lab)
       return(output_more)
    } else {
        output <- list(est_K=est_k,  OASW_est_K=oasw.value.est.k, clus_vect_est=clus_lab_est,  all_OASW = all_best_avg_silh) #, all_clus=all_best_clus_lab)
     return(output)
    }

}
 
