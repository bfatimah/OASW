#' @name pamsil
#' 
#' @title pamsil clustering
#'
#' @description An OASW clustering method based on medoids.
#'
#' 
#' @usage pamsil(dmat, distmethod = "euclidean", kmin=2, kmax=12)
#' @param dmat Either a matrix or data frame of observed values or
#' a vector of pairwise distances between observations. 
#' In first case the row represent observations to cluster and
#' columnn represents the variables. If data matrix is provided needs 
#' to specify the distance method as well. In second case usually an
#' object of class \code{"dist"}. Missing values are not allowed in both cases.
#' @param distmethod distance method to be used.
#' @param kmin minimum number of clusters for estimation of number of clusters.
#' @param kmax maximum number of clusters for estimation of number of clusters.
#'
#' @return Returns a list having following components:
#' \describe{
#' \item{est_K}{ number of clusters estimated by pamsil algorithms. All other results are based on this estimated number.}
#' \item{clus_lab}{pamsil clustering labels against estimated \code{K}.}
#' \item{silh}{optimum ASW value for each data point in the clustering.}
#' \item{clus_size}{number of observations in clusters. }
# \item{index_by_cluster}{represents the data index in each cluster. The order is same as that in \code{clus_vect}}
#' \item{avg_silh}{ASW value for pamsil clustering.}
#' \item{avg_clus_silh}{ASW for each cluster.}
#' \item{iter}{number of iteration taken by the algorithm to converge.}
#' \item{avg_silh_kmin_kmax}{ASW values for the clusterings correspondng to the number of clusters in the range \code{kmin} to \code{kmax} number of clusters.}
#'}
#' @details  This function is based on the standalone \code{C} functions written by Van et al. (2003). 
#' \code{pamsil()} sccepts both data matrix and distances.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @examples
#' dmat <- iris[,1:4]
#' dys <- dist(dmat)
#' oasw_clustering <- pamsil(dys, 2, 4)
#' oasw_clustering <- pamsil(dmat, distmethod = "manhattan", 2, 4)
#'
#'@references
#' Van der Laan, M., Pollard, K., & Bryan, J. (2003). A new partitioning around medoids
#'  algorithm. \emph{Journal of Statistical Computation and Simulation}, 73(8), 575-584.
#' @importFrom Rcpp sourceCpp evalCpp
#' @export pamsil

#---------- 
#PAMSIL 
#---------- 
pamsil <- function(dmat, distmethod = "euclidean", kmin=2, kmax=12){
  
  if(inherits(dmat, "dist") == TRUE){
    dys <- dmat
    n <- (1 + sqrt(1+8*length(dys)))/2
    dmat <- "NULL"
    indicator <- "TRUE"
  } else {
    n <- nrow(dmat)
    dys <- dist(dmat, method = distmethod)
    indicator <- "FALSE"
  }
  
  
  #estimating number of clusters
  avg_silh_k <- numeric(kmax)
  bestavg <- -2
  for(KK in kmin:kmax){
    aa_oasw <- pamsilFix(dys, KK, distmethod = "euclidean")
    avg_silh_k[KK] <- aa_oasw$avg_silh
    if(bestavg < aa_oasw$avg_silh){
      bestavg  <- aa_oasw$avg_sil
      pamsil_res <- aa_oasw
      pamsil_est_K <- KK
    }
    
  } #estimation of number of clusters ends here, so does pamsil
  
  #return(list(pamsil_res, avg_silh_k=avg_silh_k))
  #return(list(pamsil_res))
  #return(pamsil_res)
  return(list(est_K=pamsil_est_K, clus_lab = pamsil_res$clus_lab, silh = pamsil_res$silh, clus_size = pamsil_res$clus_size, avg_silh = pamsil_res$avg_silh, avg_clus_silh=pamsil_res$avg_clus_silh, iter =  pamsil_res$iter, avg_silh_kmin_kmax=avg_silh_k))
}
  
 
