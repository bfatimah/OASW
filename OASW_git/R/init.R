#' @name init
#' @title initialization function for OSIL and FOSIL. 
#'
#' @description  This function was originally wrote to provide the initialization for OSIL and FOSIL
#'  but can have standalone usage as well. More discription to follow. 
#'
#' @usage init(dmat, dys, KK=2, indicator=1)
#'@export

#################
#Initialization
#################
#indicator is true means distance between pair of observations is provided if indicator = FALSE data is provided
init <- function(dmat, dys, distmethod = "euclidean", KK, indicator) {
    
    n <- (1 + sqrt(1+8*length(dys)))/2

    Res_km <- Res_pam <- Res_avg <-  Res_ward <- Res_mod <- numeric(n)
    sil_avg <- sil_ward <- sil_pam <- sil_km <- sil_mod <- numeric(1)
    
    
    #Calculating clustering from all initialization methods from 1 - k
    #wards
    init_ward <- hclust(dys, method = "ward.D2")
    Res_ward <- cutree(init_ward, KK)
     #print("this is wards results")
     #print(Res_ward)
     #str(Res_Ward)
    sil_ward <- mean(silhouette(Res_ward, dys)[,3])
   
   #average
    init_avg <- hclust(dys, method = "average") #now this dys is normal dist() output whereas OASW takes filldys
    Res_avg <- cutree(init_avg, KK)
     #print("this is average  results")
    #print(Res_avg)
    sil_avg <- mean(silhouette(Res_avg, dys)[,3])
   
        
    #pam
    Res_pam  <- pam(dys, KK)$clustering
    sil_pam <- mean(silhouette(Res_pam, dys)[,3])

    
    #mod & kmeans #Thee is no need of showing this error if you will OSIL is not showing output for the case when in input is distance not data 
    #if(indicator == "TRUE"){
    #    skip("skipping model-based clustering while in initialization stage") #check will this message appeare on screen when you call this
    #}
    
    if(indicator == "FALSE"){
    #Kmeans
    Res_km  <- kmeans(dmat, KK)$cluster
    sil_km <- mean(silhouette(Res_km, dys)[,3])
        
    #mod
    init.res.model  <-  Mclust(dmat, KK, verbose = FALSE)
    Res_mod <- matrix(init.res.model$classification, nrow=n)   
    sil_mod <- mean(silhouette(Res_mod, dys)[,3])
    
    }

#make use of ifelse for the following aafter check.
   if(indicator == "FALSE"){
    sil_all <- data.frame(sil_avg, sil_ward, sil_km, sil_pam, sil_mod)
    sil.max <- which.is.max(sil_all)
    asw_best <-  sil_all[sil.max]
        if(sil.max == 1){
            lab_best <- Res_avg
        } else if (sil.max == 2){
            lab_best <- Res_ward
        } else if(sil.max == 3){
            lab_best <- Res_km
        }  else if (sil.max == 4){
            lab_best <- Res_pam
        }  else if (sil.max == 5){
            lab_best <- Res_mod
        }
   }
        
    
    if(indicator == "TRUE"){
            sil_all <- data.frame(sil_avg, sil_ward, sil_pam)
            sil.max <- which.is.max(sil_all)
            asw_best <-  sil_all[sil.max]
            if(sil.max == 1){
                lab_best <- Res_avg
            } else if (sil.max == 2){
                lab_best <- Res_ward
            }   else if (sil.max == 3){
                lab_best <- Res_pam
            }
        }
   
   
    if(indicator == "FALSE"){
            methodsName <- c("average", "Ward's", "kmenas", "pam", "model-based")
    } else { methodsName <- c("average", "Ward's", "pam") }
    
        
        out_init <- list(lab_best = lab_best, asw_best = asw_best, best_init_method = methodsName[sil.max])
return(out_init)

}



#Make FOSIL available for distances only. 



 
