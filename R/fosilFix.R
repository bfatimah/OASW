#' @name fosilFix
#' @title Fast osil-fixed number of clusters
#'
#' @description  This is the fast version of the OSil algorithm for a single number of cluster. 
#' OSil is an optimum average silhouette width (OASW) clustering method that 
#' donot use any kind of cluster centriods. Only data is needed as input. 
#' The algorithm can estimate number of clusters.
#'
#' @usage  fosilFix(dmat, distmethod="euclidean", k)
#'
#' @param dmat  Either a numeric matrix or data frame of observed values.
#'  The row represent observations to cluster and columnn represents the variables.
#' @param distmethod  the distance method to be used. Current available methods are "euclidean", 
#' "maximum", "manhattan", "canberra", "binary" or "minkowski". See \code{\link[stats:dist]{dist}} for more details on these methods.
#' @param k number of clusters
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
#' k <- 4
#' oasw_clus <- fosilFix(dmat, distmethod="euclidean", k)
#' plot(dmat, col = oasw_clus$clus_lab, pch = 16, cex = 1.5)
#' @references
#' Batool F. (2019). Optimum average silhouette width clustering. 
#' \emph{PhD Thesis}, University College London.
#' 
#' @export
fosilFix <- function(dmat, distmethod="euclidean", k){
  n <- nrow(dmat)
  v <- 0.2
  d <- n*v
  many.sample <- 2
  bestavgsam <- -2
  index <- 1:n
  #random sampling many.sample times
    best_asw <- -2
    for(h in 1:many.sample){       
      sam.index <- sample(index, d)
      sam.data <- dmat[sam.index,]
      #sam.dys <- dist(sam.data)
      
      out_init <- init(sam.data, k, distmethod = "euclidean")
      
      if(best_asw < out_init$asw_best){
        best_asw <- out_init$asw_best
        lab_best_sam  <- out_init$lab_best
        sam_data_best <- sam.data
        sam_index_best <- sam.index
      }
    } #sampling process ends here
    
    #clustering sample data now
    osil_sam <- osilFix(dist(sam_data_best), d, k,  lab_best_sam)
    
    colect.lab.est <- numeric(n)
    save.out.est <- numeric(k)
    
    #skipping the sample data while assigning the remainning data to clusters
    rem.index <- index[-sam_index_best]
    r <- 1
    for(i in 1:n){
      if(r <= (n-d)){
        pass.data.est  <- rbind(sam_data_best, dmat[rem.index[r], ])
      }
      
      for(h in 1:k){
        pass.lab.est <- c(osil_sam$clus_lab, h)
        dys <- dist(pass.data.est) #calculate the distance of this point only
        save.out.est[h] <- mean(silhouette(pass.lab.est, dys)[, 3])
      }
      colect.lab.est[rem.index[r]] <- which.max(save.out.est)
      r <- r + 1
    } #n for ends here. all data is assigned to clusters
    
    #filling the sample labels back in now to get the entire clustering label vector correspounding to the new data order now.
    #this will give labels according to new ordering of data
    r <- 1
    for(i in 1:n){  #to first fill the labels for the  sample data
      if(r <= d){
        colect.lab.est[sam_index_best[r]] <-  osil_sam$clus_lab[r]
        r <- r+1
      }	#outer if 
    }	#for ended
    
    #plot(data, col = 	colect.lab.est)
    
    #for k  labels are colect.lab.est, dataset and OASW value
  
  #writing output 
  #caluclate cluster_sizes
  clus_size <- table(colect.lab.est)
  
  #SW for all data points
  silh=silhouette(colect.lab.est, dist(dmat))[, 3]
  
  os_data <- mean(silhouette(colect.lab.est, dist(dmat))[, 3])
  
  #calculating avg_clus_silh
  avg_clus_silh <- numeric(k)  
  for(j in 1:k){
    for(i in 1:n){
      if(colect.lab.est[i] == j){
        avg_clus_silh[j] <- silh[i]
      }
    }
  }
  
  out <- list(n=n, k=k,  clus_lab=colect.lab.est, clus_size=as.integer(clus_size), silh=silh, avg_clus_silh=avg_clus_silh, avg_silh=os_data)
  
  return(out)
} #fosilFix function ends here.  



