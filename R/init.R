#' @name init
#' @title initialization function for OASW clustering algorithms
#'
#' @description  This function was originally wrote to provide the initialization for osil
#'  and fosil but can have standalone usage.  The function takes a fixed number of clusters.
#'
#' @usage init(dmat, K, distmethod = "euclidean", ...)
#' @param dmat either a numeric matrix or data frame of observed values or 
#' pairwise distances between observations. 
#' @param K number of cluster
#' @param distmethod distance method to be used for the calculation of
#'  pairwise distances between observations.
#' @param ... additional parameters
#' @return Returns a list having following components
#' \describe{
#' \item{n}{number of observations.}
#' \item{K}{number of clusters.}
#' \item{lab_best}{clustering labels corresponding to the best ASW clustering among the  initialization methods. }
#' \item{asw_best }{best ASW value among the initialization methods. }
#' \item{best_init_method}{name of the best initialization methods based on ASW value. }
#' }
#' @details  The function is originally written as an initialization function for osil and fosil 
#' clusterings. The \code{init} functions returns the best clustering out of the six 
#' methods initialization methods based on the ASW values. Several clustering methods were considered
#' in a systematic simulation set-up to find out the best for the OASW clustering initialization.
#'  Among those considered six showed good performance for single or multiple data generating structures 
#'  considered in the simulations (see Batool, 2019 for results and discussion).
#' The six best initialization methods are average (Sokal and Michener 1958), Ward's  (Ward 1963), pam (Kaufman and Rousseeuw 1990), 
#' kmeans (Hartigan and Wong 1979), and model-based (Fraley and Raftery 1998) clustering. Works with both distances and data matrix.
#' Currently, in case distances are provided \code{kmeans} and \code{Mclust} are
#' excluded from initialization methods.
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#' 
#' @seealso \code{\link{hclust}}, \code{\link{pam}}, \code{\link{kmeans}}, \code{\link{Mclust}}
#' functions to pass additional argumnets to \code{\link{init}}.
#' 
#' @examples
#' dmat <- TwoGaussian(100)
#' dys <- dist(dmat$data)
#' plot(dmat$data, col = dmat$truelab)
#' init_res_1 <- init(dmat$data, KK=2, distmethod = "euclidean")
#' init_res_2 <- init(dys, KK=2, distmethod = "euclidean")
#' print(init_res_1)
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
#' @import mclust
#' @importFrom nnet which.is.max
#' @export

#-----------------
#Initialization
#-----------------
#indicator=TRUE: distances between pairs of observations are provided,
#if indicator=FALSE data matrix is provided
init <- function(dmat, K, distmethod = "euclidean", ...) {
    
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

    Res_km <- Res_pam <- Res_avg <- Res_single <- Res_ward <- Res_mod <- numeric(n)
    sil_avg <- sil_ward <- sil_single <- sil_pam <- sil_km <- sil_mod <- numeric(1)
    
    #Calculating clustering from all initialization methods from 1-k
    #wards
    init_ward <- hclust(dys, method = "ward.D2")
    Res_ward <- cutree(init_ward, K)
    sil_ward <- mean(silhouette(Res_ward, dys)[,3])
   
    #average
    init_avg <- hclust(dys, method = "average")
    Res_avg <- cutree(init_avg, K)
    sil_avg <- mean(silhouette(Res_avg, dys)[,3])
    
    #single
    init_single <- hclust(dys, method = "single") 
    Res_single <- cutree(init_single, K)
    sil_single <- mean(silhouette(Res_single, dys)[,3])
           
    #pam
    Res_pam  <- pam(dys, K)$clustering
    sil_pam <- mean(silhouette(Res_pam, dys)[,3])

    if(indicator == "FALSE"){
    #kmeans
    Res_km  <- kmeans(dmat, K)$cluster
    sil_km <- mean(silhouette(Res_km, dys)[,3])
     
    #model-based
    init.res.model  <-  Mclust(dmat, K, verbose = FALSE)
    Res_mod <- matrix(init.res.model$classification, nrow=n)   
    sil_mod <- mean(silhouette(Res_mod, dys)[,3])
    
    }

   if(indicator == "FALSE"){
    sil_all <- cbind(sil_avg, sil_ward, sil_single, sil_km, sil_pam, sil_mod)
    sil.max <- which.is.max(sil_all)
     asw_best <- numeric(0)
    asw_best <-  sil_all[sil.max]
        if(sil.max == 1){
            lab_best <- Res_avg
        } else if (sil.max == 2){
            lab_best <- Res_ward
        } else if(sil.max == 3){
            lab_best <- Res_single
        }  else if (sil.max == 4){
            lab_best <- Res_km
        }  else if (sil.max == 5){
            lab_best <- Res_pam
        }  else if (sil.max == 6){
            lab_best <- Res_mod
        }
   }
        
    
    if(indicator == "TRUE"){
            sil_all <- data.frame(sil_avg, sil_ward, sil_single, sil_pam)
            sil.max <- which.is.max(sil_all)
            asw_best <-  sil_all[sil.max]
            if(sil.max == 1){
                lab_best <- Res_avg
            } else if (sil.max == 2){
                lab_best <- Res_ward
            }   else if (sil.max == 3){
                lab_best <- Res_single
            }   else if (sil.max == 4){
                lab_best <- Res_pam
            }
        }
   
    if(indicator == "FALSE"){
            methodsName <- c("average", "Ward's", "Single", "kmenas", "pam", "model-based")
    } else { methodsName <- c("average", "Ward's", "Single", "pam") }
    
        
    out_init <- list(n=n, K=K, lab_best = lab_best, asw_best = asw_best, best_init_method = methodsName[sil.max])
    return(out_init)

}


 
