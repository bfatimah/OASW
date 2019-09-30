#' @title OASW_fun
#' @description This is a wrapper function for C functions.  This function is called from within OSIL funciton.  Can be called standalone for clustering for fixed number of clusters with an intialization label set. This will return an OASW clsutering based on labels
#'
#' @param dys pairwise distances between observations
#' @param n number of observations
#' @param K number of clusters
#' @param clus_lab initial clustering labels
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
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}


################################
# "OASW2" baesd on labels.
################################

OASW_fun <- function(dys, n, K, clus_lab){
    #n is no. of observations
    #K is no. of clusters
    #Initializing algorithm
    
    #METHODS <- c("Rand", "SRSO", "ward.D", "single", "complete", "average", "mcquitty", "ward.D2", "PAM", "k-means", "Model", "Spectral", "known")
    #check_aval_meth <- pmatch(init.method, METHODS)
    # if (is.na(check_aval_meth))
    #   stop("invalid clustering method", paste("", init.method))
    #if (check_aval_meth == -1)
    #  stop("ambiguous clustering method", paste("", init.method))
    
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
    
    output <- list(n = n, K = K, clus_lab = alt_clus_lab, clus_size = clus_size, silh = silh, avg_clus_silh = avg_clus_silh, avg_silh = avg_silh, iter = iter )
    output
}

