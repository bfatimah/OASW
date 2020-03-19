#' @name osilFix
#' @title osil clustering-fixed number of clusters
#' @description Produces a clustering solution for a fixed number of clusters
#' @encoding UTF-8
#' @param dys pairwise distances between observations.
#' @param n number of observations.
#' @param K number of clusters.
#' @param clus_lab clustering labels.
#' @return Returns a list having following components.
#' \describe{
#' \item{n}{total number of data points.}
#' \item{K}{estimated number of clusters.}
#' \item{clus_lab}{clustering labels.}
#' \item{clus_size}{cluster sizes.}
#' \item{silh}{ASW value of each cluster.}
#' \item{avg_clus_silh}{ASW for each cluster.}
#' \item{avg_silh}{ASW value for the clustering.}
#' \item{iter}{number of iteration taken by the algorithm to converge.} }
#' 
#' @details This is a wrapper function for C++ functions. 
#' This function is called from within osil function.  
#' Can be called standalone for clustering for fixed 
#' number of clusters with an intialization clustering label set. 
#' This will return an OASW clustering based on labels.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#' @importFrom Rcpp sourceCpp evalCpp
#' @export


################################
# "OASW2" baesd on labels.
# Used for osil and fosil and fosilFix
################################
#dys is pairwise distance between observations
osilFix <- function(dys, n, K, clus_lab){
    #n is no. of observations
    #K is no. of clusters
    #Initializing algorithm
    
    clus_lab <- clus_lab - 1
    alt_clus_lab <- integer(n)
    clus_size <- integer(K)
    iter <- integer(1)
    dys_i <- double(n)
    avg_dys_clus <- double(n*K)
    silh <- double(n)
    altsilh <- double(n)
    avg_clus_silh <- double(K)
    avg_clus_silhtwo <- double(K)
    avg_silh <- double(1)
    disty <- numeric(n*(n-1)/2)
    disty <- filldys(dys)
    
    sil_lab_swap(K, n, clus_lab, alt_clus_lab, clus_size, disty, iter, dys_i, avg_dys_clus, silh, altsilh)
    
    clustyanlys(K, n, alt_clus_lab, clus_size, silh, avg_clus_silh, avg_clus_silhtwo, avg_silh)
    
    clus_lab <- clus_lab + 1
    alt_clus_lab <- alt_clus_lab + 1
    
    output <- list(n = n, K = K, clus_lab = alt_clus_lab, clus_size = clus_size, silh = silh, avg_clus_silh = avg_clus_silh, avg_silh = avg_silh, iter = iter)
    output
}

