#' @name pairplots 
#' @title  pairplots 
#' @description pair plots for more than two demensional data
#' 
#' @usage pairplots(data, labels)
#'
#' @param data data set to plot.
#' @param labels labels against which plotting is needed.
#' @return Returns a plot.
#' @details none
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- C7D10(70)
#' pairplots(dmat$data, dmat$truelab)   
#  @import GGally
#' @export


pairplots <- function(data, labels){
  data <- data.frame(data)
  truelab <- factor(labels)
  model <- ggpairs(data,  aes(colour = truelab, alpha = 0.5),  
                   upper = list(continuous = "density"),
                   lower = list(continuous = "points"), 
                   diag = list(continuous =  "barDiag"),  
                   xlab = NULL, ylab = NULL, axisLabels = c("show"),
                   cardinality_threshold = 15) + 
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
    )  
  return(model)
}

 

