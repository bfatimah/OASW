#' @name FOSIL
#' @title Fast optimum avearage silhouette width clustering
#'
#' @description  This is the fast version of the OSil algorithm. 
#' OSil is an optimum average silhouette width (OASW) clustering method that 
#' donot use any kind of cluster centriods. Only data is needed as input. 
#' The algorithm can estimate number of clusters.
#'
#' @usage FOSIL(dmat, distmethod="euclidean", k=4, v=0.2, many.sample=2)
#'
#' @param dmat  Either a numeric matrix or data frame of observed values.
#'  The row represent observations to cluster and columnn represents the variables.
#' @param distmethod  the distance mehtod to be used. Current available methods are "euclidean", 
#' "maximum", "manhattan", "canberra", "binary" or "minkowski". 
#' See \code{\link[stats:dist]{dist}} for more details on these methods.
#' @param k maximum value to be used for the estimation of number of clusters.
#' @param v data proportion to be used for samplig.
#' @param many.sample how many times random samples of size v to be taken from data.
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
#' \item{iter}{number of iteration taken by the algorithm to converge.}
#' \item{best_init_method}{method that gave the best initialization i.e., maximum ASW.}
#' }
#' @details  OSIL has an initialization phase. Based on extensive simulaitons the best
#'  initialization methods from among a wide range of existing clustering methods, 
#'   for the algorithm has been identified namely average, Ward's, pam,
#'    kmeans, and model-based clustering. In case distances are provided for 
#'    clustering kmeans and model-based clustering are excluded from the initialization methods.
#'
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#' @seealso \code{\link{OSIL}} for not so fast version.
#'
#' @examples
#' require(mlbench)
#' dmat <- mlbench.shapes(100)$x
#' oasw_clus <- FOSIL(dmat) 
#' plot(dmat, col = oasw_clus$clus_lab, pch = 16, cex = 1.5)
#'@references
#' Hartigan, J. A. and Wong, M. A. (1979). Algorithm AS 136: A K-means clustering algorithm. 
#' \emph{Applied Statistics}, 28, 100–108.
#' 
#' Kaufman, L. and P. J. Rousseeuw (1990). Finding groups in data: 
#' an introduction to cluster analysis, Volume 344. John Wiley and Sons.
#' 
#' R. R. Sokal and C. D. Michener. A statistical method for evaluating systematic relationships.
#'  \emph{Univesity Kansas Science Bulletin}, 38(22):1409–1438, 1958.
#'  
#' J. H. Ward Jr. Hierarchical grouping to optimize an objective function.
#'  \emph{Journal of the American Statistical Association}, 58(301):236–244, 1963.
#'  
#' C. Fraley and A. E. Raftery. How many clusters? which clustering method?
#'  answers via model-based cluster analysis. \emph{The Computer Journal}, 41(8):578–588, 1998.
#'  
#' Batool F. (2019). Optimum average silhouette width clustering. 
#'  \emph{PhD Thesis}, University College London.
#' 
#' @export
FOSIL <- function(dmat, distmethod="euclidean", k=12, v=0.2, many.sample=2){
n <- nrow(dmat)
d <- n*v
bestavgsam <- -2
Res_lab <- matrix(0, nrow=n, ncol=k)
os_data <- numeric(k)
index <- 1:n
for(KK in 2:k){
    best_asw <- -2
    for(h in 1:many.sample){       
       sam.index <- sample(index, d)
       sam.data <- data[sam.index,]
       sam.dys <- dist(sam.data)
       
       indicator = "FALSE"  
       out_init <- init(sam.data, sam.dys, distmethod = "euclidean", KK, indicator)
                
        if(best_asw < out_init$asw_best){
            best_asw <- out_init$asw_best
            lab_best_sam  <- out_init$lab_best
            sam_data_best <- sam.data 
            sam_index_best <- sam.index
        }
    } #sampling process ends here
    #plot(sam_data_best, col =  lab_best_sam)
    #clustering sample data now
    osil_sam <- OASW_fun(dist(sam_data_best), d, KK,  lab_best_sam)
 #plot(sam_data_best, col =  osil_sam$clus_lab)
 colect.lab.est <- numeric(n)
 save.out.est <- numeric(KK)
 
 #skipping the sample data while assigning the remainning data to clusters
rem.index <- index[-sam_index_best]
r <- 1
for(i in 1:n){
	if(r <= (n-d)){
	pass.data.est  <- rbind(sam_data_best, data[rem.index[r], ])
    }
    
    for(h in 1:KK){
            pass.lab.est <- c(osil_sam$clus_lab, h)
            dys <- dist(pass.data.est) #calculate the distance of this point only
            save.out.est[h] <- mean(silhouette(pass.lab.est, dys)[, 3])
        }
colect.lab.est[rem.index[r]] <- which.max(save.out.est)
    r <- r + 1
} #n for ends here

#filling the sample labels backin now to get the entire clustering vector in the original data order
    r <- 1
	for(i in 1:n){  #to first fill the labels for the 
	if(r <= d){
		colect.lab.est[sam_index_best[r]] <-  osil_sam$clus_lab[r]
		r <- r+1
	}	#outer if 
	}	#for ended
	
	#plot(data, col = 	colect.lab.est)
	
#for KK store labels, dataset and OASW value
    for(i in 1:n){
        Res_lab[,KK] <- matrix(colect.lab.est, nrow = n)
    }
    os_data[KK] <- mean(silhouette(Res_lab[,KK], dist(data))[, 3])   
} #KK for ended 
   
  #estimation of K
osj_K <- which.max(os_data)
 Res_lab[,osj_K]
 #plot(data, col =  Res_lab[,osj_K])
  
  #writing output 
  #caluclate cluster_sizes
  clus_size <- table(Res_lab[,osj_K])
  
  #SW for all data points
  silh=silhouette(Res_lab[,osj_K], dist(data))[, 3]
  
  #calculating avg_clus_silh
  avg_clus_silh <- numeric(osj_K)
 for(j in 1:osj_K){
 	for(i in 1:n){
 		if(Res_lab[i,osj_K] == j){
 			avg_clus_silh[j] <- silh[i]
 		}
 	}
 }
 
  FOSIL_OUT <- list(n=n, K=osj_K,  clus_lab=Res_lab[,osj_K], clus_size= clus_size, silh=silh, avg_clus_silh=avg_clus_silh, avg_silh=os_data[osj_K])
 
  return(FOSIL_OUT)
} #FOSIL function ends here.  
	
	
	
