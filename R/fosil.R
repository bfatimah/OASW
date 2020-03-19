#' @name fosil
#' @title Fast osil-estimation of number of clusters
#'
#' @description  This is the fast version of the OSil algorithm. 
#' OSil is an optimum average silhouette width (OASW) clustering method that 
#' donot make use of any kind of cluster centriods. Only data is needed as input. 
#' The algorithm can estimate number of clusters.
#'
#' @usage fosil(dmat, distmethod="euclidean", kmin=2, kmax=12)  
#'
#' @param dmat  Either a numeric matrix or data frame of observed values.
#'  The row represent observations to cluster and columnn represents the variables.
#' @param distmethod  the distance method to be used. Current available methods are "euclidean", 
#' "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' See \code{\link[stats:dist]{dist}} for more details on these methods.
#' @param kmin minimum number of clusters
#' @param kmax maximum value to be used for the estimation of number of clusters.
# @param v data proportion to be used for samplig.
# @param many.sample how many times random samples of size v to be taken from data.
#'
#' @return Returns a list having following components:
#' \describe{
#' \item{n}{total number of data points.}
#' \item{K}{estimated number of clusters.}
#' \item{clus_lab}{clustering label vector.}
#' \item{clus_size}{cluster sizes.}
#' \item{silh}{OASW value of each cluster.}
#' \item{avg_clus_silh}{OASW for each cluster.}
#' \item{avg_silh}{OASW value for the clustering.}
#' \item{$avg_silh_kmin_kmax}{ASW value for all the clusterings.}
#' 
# \item{iter}{number of iteration taken by the algorithm to converge.}
# \item{best_init_method}{method that gave the best initialization i.e., maximum ASW.}
#' }
#' @details  osil has an initialization phase. Based on extensive simulaitons the best
#'  initialization methods from among a wide range of existing clustering methods, 
#'   for the algorithm has been identified namely average, Ward's, pam,
#'    kmeans, and model-based clustering. In case distances are provided for 
#'    clustering kmeans and model-based clustering are excluded from the initialization methods.
#'
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @seealso \code{\link{osil}} for not so fast version.
#'
#' @examples
#' require(mlbench)
#' dmat <- mlbench.shapes(100)$x
#' oasw_clus <- fosil(dmat) 
#' oasw_clus <- fosil(dmat, distmethod="euclidean", kmin=2, kmax=5)
#' plot(dmat, col = oasw_clus$clus_lab, pch = 16, cex = 1.5)
#' @references
#' Batool F. (2019). Optimum average silhouette width clustering. 
#' \emph{PhD Thesis}, University College London.
#' 
#' @export
fosil <- function(dmat, distmethod="euclidean", kmin=2, kmax=12){
  n <- nrow(dmat)
  v <- 0.2
  d <- n*v
  avg_silh_k <- numeric(kmax)
  #estimation of k
  bestavg <- -2
  for(KK in kmin:kmax){
    out_fosilFix <- fosilFix(dmat, distmethod="euclidean", KK)
    avg_silh_k[KK] <- out_fosilFix$avg_silh
    if(bestavg < out_fosilFix$avg_silh){
      bestavg  <- out_fosilFix$avg_sil
      best_res <- out_fosilFix
    }
  } #KK for ended 
  
  out <- list(n=best_res$n, est_k=best_res$k,  clus_lab=best_res$clus_lab, clus_size= best_res$clus_size, silh=best_res$silh, avg_clus_silh=best_res$avg_clus_silh, avg_silh=best_res$avg_silh, avg_silh_kmin_kmax= avg_silh_k) #iter=best_res$iter
  
  return(out)
} #FOSIL function ends here.  


