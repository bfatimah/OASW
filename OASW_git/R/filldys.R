#' @title filldys
#'
#' @description An auxilary function for HOSIL, OSIL, FOSIL and PAMSIL.
#'
#' @param dys A vector of pairwise distances between observations. Usually an object of class "dist"
#' @return Returns a vector of pairwise distances between observations. The total entries are \code{n*(n-)/2+1}
#' @details filldys returns a vector of distances that append 0 in the start. This is to deal with the main diagonal entry of the distance matrix. 0 is appended in the begining to obtain the diagonal values which are all zeros in the distance matrix whenever needed. Note that other alternatives might not work for instance c(0, dys). Do not try to change


################################
## filldys
################################
filldys <- function(dys){ #mdc = as.matrix(dist) only
    n <- (1 + sqrt(1+8*length(dys)))/2
    mdc <- as.matrix(dys)
    n_dys <- n*(n-1)/2
    vect <- double(n_dys)
    count = 1
    
    for(i in 1:n){ #this is second index
        for(j in 1:n){ #this is first index
            if(j < i){
                vect[count] <-  mdc[i, j]
                count = count + 1
            }
        }
    }
    
    n_dysnn <- n*(n-1)/2+1
    dys <- double(n_dysnn)
    for(i in 1:n_dysnn)
    {
        if(i == 1){
            dys[i] = 0
        } else {
            dys[i] = vect[i-1]
        }
    }
    return(dys)
}

