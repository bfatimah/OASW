#' @name plot2d
#' @title Enhanced 2d plotting
#' @description plots for two demensional data
#' 
#' @usage plot2d(data, labels)
#'
#' @param data data to plot
#' @param labels labels against which plotting is needed
#' @return Returns a plot
#' @details none
#' @author Fatima Batool \email{ucakfba@ucl.ac.uk}
#'
#'
#' @examples
#' dmat <- C7D10(1000)
#' plot2d(dmat$data, dmat$truelab) 
#  @import ggplot2
#' @export


plot2d <- function(data, labels){
  data <- data.frame(data)
  model <- ggplot(data, aes(x = data[,1], y = data[,2])) + 
    geom_point(color =  labels, size = 2) + 
    theme_bw() +labs(x ="", y = "") + theme(text = element_text(size=12))  + 
    theme(panel.border = element_rect(colour = "black", size =1.5), 
          axis.title = element_text(face = "bold",size = rel(1)),	
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(size = 12),
          axis.ticks = element_line(),
          panel.grid.major = element_blank(), #element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          plot.margin=unit(c(10,5,5,5),"mm"))  
  return(model)
}

