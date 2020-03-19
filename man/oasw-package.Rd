\name{oasw-package}
\alias{oasw-package}
\docType{package}
\title{
  Optimum average silhouette width clustering methods
}
\description{
  The package implements the hierarchical and partitional clustering methods
  based on the optimization of the average silhouette width index.  The following are the major functions of the package.
}
\section{Clustering algorithms}{
 \describe{
 \item{pamsil}{Computes PAMSIL clustering as introduced
    in Van der Laan (2003).}
 \item{osil}{Computes osil clustering as introduced
    in Batool (2019).}
 \item{fosil}{Computes fosil clustering as introduced
    in Batool (2019).}
  \item{hosil}{Computes hosil clustering as introduced
    in Batool (2019).}
 }
 }
\author{
Fatima Batool

Maintainer: Fatima Batool <fatima.batool.14@ucl.ac.uk>
}
\references{
  Batool F. (2019). Optimum average silhouette width clustering.
  \emph{PhD Thesis}, University College London.
  
  Batool, F., (2019). A new hierarchical clustering algorithm based on 
  optimization of ASW linkage criterion. \emph{Journal of classification}.
  \url{https://arxiv.org/abs/1909.12356}.
  
  Batool, F., and Hennig, C. (2019). Characterization and Development of 
  Average Silhouette Width Clustering. \emph{Advanced in data analysis and classification}.
  \url{https://arxiv.org/abs/1910.11339}.  
  
  Batool, F., and Hennig, C. (2019). Initializations and related challenges for clustering 
  by optimizing the average silhouette width.  \url{https://arxiv.org/abs/1910.08644}.    
 
 Van der Laan, M., Pollard, K., & Bryan, J. (2003). A new partitioning around medoids
  algorithm. \emph{Journal of Statistical Computation and Simulation}, 73(8), 575-584.
}
\keyword{ Rcpp }
\seealso{
  \code{\link{kmeans}} \code{\link{pam}} \code{\link{hclust}} \code{\link{Mclust}}
}
