#' @name OSIL
#' 
#' @title Optimum avearge silhouette width clustering
#'
#' @description  An OASW clustering method that does not use any kind of cluster centriods.
#'
#' 
#' @usage OSIL(dmat, distmethod = "euclidean", k=4)
#' 
#' @param dmat  either a numeric matrix or data frame of observed values or pairwise distances between observations. In first case the row represent observations to cluster and columnn represents the variables. If data matrix is provided needs to specify the distance method as well. In second case usually an object of class \code{dist}. Missing values are not allowed in both cases.
#' @param distmethod  the distance mehtod to be used. Current available methods are "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". See \code{\link[stats:dist]{dist}} for more details on these methods.
#' @param k maximum value to be used for the estimation of number of clusters
#'
#' @return Returns a list having following components
#' \describe{
#' \item{n}{total number of data points}
#' \item{K}{estimated number of clusters}
#' \item{clus_lab}{clustering label vector}
#' \item{clus_size}{cluster sizes}
#' \item{silh}{OASW value of each cluster }
#' \item{avg_clus_silh}{OASW for each cluster}
#' \item{avg_silh}{OASW value for the clustering}
#' \item{iter}{number of iteration taken by the algorithm to converge} }
#' @details  The data given in matrix form is clustered by the newely proposed OSil algorithm based on the optimization of ASW index. OSIL has an initialization phase. Based on extensive simulaitons the best initialization methods from among a wide range of existing clustering methods,  for the algorithm has been identified namely average ref, Ward's  ref, pam ref, kmeans ref and model-based ref clustering. In case distances are provided for clustering kmeans  and model-based clustering are excluded from the initialization methods.
#'
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @examples
#' require(mlbench)
#' dmat <- mlbench.shapes(100)$x
#' dys <- dist(dmat)
#' oasw_clus <- OSIL(dmat) 
#' oasw_clus <- OSIL(dys)
#' plot(dmat, col = oasw_clus$clus_lab, pch = 16, cex = 1.5)
#'@references
#' Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: 
#' A K-means clustering algorithm. \emph{Applied Statistics}, 28, 100–108.
#' 
#' Kaufman, L. and P. J. Rousseeuw (1990). Finding groups in data: 
#' an introduction to cluster analysis, Volume 344. John Wiley & Sons.
#' 
#' R. R. Sokal and C. D. Michener. A statistical method for evaluating systematic 
#' relationships. \emph{Univesity Kansas Science Bulletin}, 38(22):1409–1438, 1958.
#' 
#' J. H. Ward Jr. Hierarchical grouping to optimize an objective function. 
#' \emph{Journal of the American Statistical Association}, 58(301):236–244, 1963.
#' 
#' C. Fraley and A. E. Raftery. How many clusters? which clustering method? 
#' answers via model-based cluster analysis. \emph{The Computer Journal}, 41(8):578–588, 1998.
#' 
#' Batool F. (2019). Optimum average silhouette width clustering.  
#' \emph{PhD Thesis}, University College London.
#'  
#' @export
####################################
#OSil_1 #from Final_algorithms file
####################################
#dmat is the  data matrix of size n*p. rows are observations to clusters and columns are features or variables
#number of clusters are estimated from the range 1:k=12
OSIL <- function(dmat, distmethod = "euclidean", k=12){
 
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
 
#p <- ncol(dmat) #number of dimensions/variables in data. not needed 


#estimating number of clusters
    bestavg <- -2
    for(KK in 2:k){
        out_init <- init(dmat, dys, distmethod = "euclidean", KK, indicator)
        aa_oasw <- OASW_fun(dys, n, KK,  out_init$lab_best)
    
        if(bestavg < aa_oasw$avg_silh){
            bestavg  <- aa_oasw$avg_sil
            osilpre_res <- aa_oasw
        }
        
    } #estimation of number of clusters ends here, so does OSIL
            
return(osilpre_res)
}
 

#library(devtools)
#library(roxygen2)
#alpha <- as.package("/Users/fatimabatool/Documents/New/OASW")
#load_all(alpha)  
#document(alpha)
#devtools::check(manual=TRUE)

