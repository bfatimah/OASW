#' @name pamsilFix
#' @title  A new partitioning around medoids algorithm
#'
#' @description  An OASW clustering method based on medoids. 
#' A PAM like clustering algorithm for the optimization of ASW based on medoids proposed in Van et al. (2003).he
#'
#' @param dys  Either a matrix or data frame of observed values or a vector of pairwise distances between observations. In first case the row represent observations to cluster and columnn represents the variables. If data matrix is provided needs to specify the distance method as well. In second case usually an object of class \code{"dist"}. Missing values are not allowed in both cases.
#' @param K number of clusters
#' @param distmethod distance method to be used 
#'
#'
#' @return Returns a list having following components:
#' \describe{
#' \item{clus_lab}{pamsil clustering labels.}
#' \item{silh}{ASW value for each data point in the clustering.}
#' \item{clus_size}{number of observations in clusters.}
# \item{index_by_cluster}{represents the data index in each cluster. The order is same as that in \code{clus_vect}.}
#' \item{avg_clus_silh}{ASW for each cluster.}
#' \item{avg_silh}{ASW value for the clustering.}
#' \item{iter}{number of iteration taken by the algorithm to converge.}
#'}
#' @details  This function is based on the standalone \code{C} functions written by Van et al. (2003). Use of this 
#' function is recommended if number of clusters are known. \code{pamsilFix} accepts both data matrix and distances.
#'
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @examples
#' dmat <- iris[,1:4]
#' dys <- dist(dmat)
#' oasw_clustering <- pamsilFix(dys, 3)
#' oasw_clustering <- pamsilFix(dmat, 3, distmethod = "manhattan")
#'
#'@references
#' Van der Laan, M., Pollard, K., & Bryan, J. (2003). A new partitioning around medoids
#'  algorithm. \emph{Journal of Statistical Computation and Simulation}, 73(8), 575-584.
#' @importFrom Rcpp sourceCpp evalCpp
#' @export
#---------
# pamsil
#---------

pamsilFix <- function(dys, K, distmethod = "euclidean"){
    
    if(inherits(dys, "dist") == FALSE){
        dys = dist(dys, method = distmethod)
    }

    n <- (1 + sqrt(1+8*length(dys)))/2
    rep_ind <- integer(n)
    dys_1 <- double(n)
    dys_2 <- double(n)
    effect <- double(n)
    max_dys <- 2e-308
    n_dys <- 1+n*(n-1)/2
    #dys <- numeric(n_dys) #this is not needed
    dys <-  filldys(dys)
    
    for (i in 1:n_dys)
    {
        max_dys = max(max_dys, dys[i])
    }
    for (i in 1:n)
    {
        dys_1[i] = 1.1*max_dys +1;
    }
    
    
    #bswap contains built and swap phases of PAM
    bswap(K, n, rep_ind, dys_1, dys_2, effect, dys, max_dys)
    
    #Calling silswap. This is optimizing average silhouette width (OASW). OASW is based on swapping
    medoid <- integer(K)
    altmeds <- integer(K)
    clus_vect <-  integer(n)
    clus_size <- integer(K)
    dys_i <-  double(n)
    avg_dys_clust <-  double(n*K)
    silh <-  double(n)
    altsilh <- double(n)
    iter <- integer(1)
    silswap(K, n, rep_ind, medoid, altmeds, clus_vect, clus_size, dys, dys_1, dys_i, avg_dys_clust, silh, altsilh, iter)
    
    #Calling clusanalysis
    avg_dys_clust <- double(n*K)
    clus_diam = double(K)
    clus_sep = double(K)
    clus_dys_1_avg = double(K)
    clus_dys_1_max = double(K)
    index_by_cluster = integer(n)
    a = double(n)
    b = double(n)
    avg_clus_silh = double(K)
    avg_silh = double(1)
    
    clusanal(K, n, rep_ind, dys_1, dys, max_dys, medoid, clus_size, clus_diam, clus_sep,
    clus_dys_1_avg,  clus_dys_1_max, index_by_cluster, clus_vect, avg_dys_clust, a, b, silh,  avg_clus_silh, avg_silh)
    
    output <- list(clus_lab = clus_vect, silh = silh, clus_size = clus_size,
    avg_clus_silh = avg_clus_silh, avg_silh = avg_silh, iter = iter) # index_by_cluster = index_by_cluster,
    output
}

#may return silently the original dys, number of observations etc.
#See Dist function in package amap by Antine Lucas it has more dist options

