#' @name osil
#' 
#' @title Optimum Average Silhouette Width clustering 
#'
#' @description  An OASW clustering method that does not use cluster centriods with estimation of number of clusters.
#'
#' @usage osil(dmat, distmethod = "euclidean", kmin=2, kmax=12)
#' 
#' @param dmat  either a numeric matrix or data frame of observed values or pairwise distances between observations. 
#' In first case the row represent observations to cluster and columnn represents the variables. 
#' If data matrix is provided then the distance method can be specfied as well. In second case usually an object of class \code{dist}.
#' Missing values are not allowed in both cases.
#' @param distmethod  the distance mehtod to be used. Current available methods are "euclidean", "maximum", "manhattan", "canberra", 
#' "binary" or "minkowski". See \code{\link[stats:dist]{dist}} for more details on these methods.
#' @param kmin minimum value to be used for the estimation of number of clusters.
#' @param kmax maximum value to be used for the estimation of number of clusters.
#'
#' @return Returns a list having following components;
#' \describe{
#' \item{n}{total number of data points.}
#' \item{K}{estimated number of clusters.}
#' \item{clus_lab}{clustering labels.}
#' \item{clus_size}{number of observations in clusters.}
#' \item{silh}{ASW value of each cluster.}
#' \item{avg_clus_silh}{\code{ASW} for each cluster.}
#' \item{avg_silh}{\code{ASW} value for the clustering for \code{K} clusters.}
#' \item{iter}{number of iteration taken by the algorithm to converge.}
#' \item{avg_silh_k}{ASW value for the clustering for \code{kmin} to \code{kmax} clusterings}.  }
#' @details  The data given in matrix form is clustered by the newely proposed OSil algorithm (see Batool 2019) based on the optimization of ASW index. 
#' osil has an initialization phase. Based on extensive simulaitons the best initialization methods from among a wide range of existing
#' clustering methods,  for the algorithm has been identified namely average (Sokal and Michener 1958), Ward's  (Ward 1963), pam (Kaufman and Rousseeuw 1990), 
#' kmeans (Hartigan and Wong 1979), and model-based (Fraley and Raftery 1998) clustering.
#' In case distances are provided for clustering kmeans and model-based clustering are excluded from the initialization methods.
#' @details \code{kmin} and \code{kmax} define the range for the estimation of number of clusters.
#' If a single valued clustering solution say \code{K} is required, specify \code{kmin=K} 
#' and \code{kmax=K} or use \code{\link{osilFix}}.
#'
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#' 
#' @examples
#' dmat <- TwoGaussian(100)$data
#' oasw_clustering <- osil(dmat)
#' dys <- dist(dmat)
#' oasw_clustering <- osil(dys)
#' plot(dmat, col = oasw_clus$clus_lab, pch = 16, cex = 1.5)
#' @references 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#' 
#' C. Fraley and A. E. Raftery (1998). How many clusters? which clustering method? 
#' answers via model-based cluster analysis. \emph{The Computer Journal}, 41(8):578 588.
#' 
#' Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: 
#' A K-means clustering algorithm. \emph{Applied Statistics}, 28, 100 108.
#' 
#' J. H. Ward Jr. (1963). Hierarchical grouping to optimize an objective function. 
#' \emph{Journal of the American Statistical Association}, 58(301):236 244.
#' 
#' Kaufman, L. and P. J. Rousseeuw (1990). Finding groups in data: 
#' an introduction to cluster analysis, Volume 344. John Wiley & Sons.
#' 
#' McQuitty, L. L. (1957). Elementary linkage analysis for isolating orthogonal and oblique 
#' types and typal relevancies. Educational and PsychologicalMeasurement 17(2), 207 229.
#'  
#' R. Sokal and C. D. Michener (1958). A statistical method for evaluating systematic 
#' relationships. \emph{Univesity Kansas Science Bulletin}, 38(22):1409 1438.
#' 
#' @importFrom stats dist hclust cutree kmeans
#' @importFrom utils data
#' @importFrom cluster silhouette pam
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom mclust Mclust mclustBIC
#' @importFrom nnet which.is.max
#' @export

#---------
#OSil 
#---------
osil <- function(dmat, distmethod = "euclidean", kmin=2, kmax=12){
    
    if(inherits(dmat, "dist") == TRUE){
        dys <- dmat
        n <- (1 + sqrt(1+8*length(dys)))/2
    } else {
        n <- nrow(dmat)
        dys <- dist(dmat, method = distmethod)
    }
    
    avg_silh_k <- numeric(kmax)
    bestavg <- -2
    for(KK in kmin:kmax){
        out_init <- init(dmat, KK, distmethod = "euclidean")
        aa_oasw <- osilFix(dys, n, KK,  out_init$lab_best)
        avg_silh_k[KK] <- aa_oasw$avg_silh
        if(bestavg < aa_oasw$avg_silh){
            bestavg  <- aa_oasw$avg_sil
            osilpre_res <- aa_oasw
        }
        
    } 
            
return(list(n = osilpre_res$n, est_K = osilpre_res$K, clus_lab = osilpre_res$clus_lab, clus_size = osilpre_res$clus_size, silh = osilpre_res$silh, avg_clus_silh = osilpre_res$avg_clus_silh, avg_silh=osilpre_res$avg_silh, iter = osilpre_res$iter, avg_silh_kmin_kmax=avg_silh_k))
}

