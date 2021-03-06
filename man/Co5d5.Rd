% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Co5d5.R
\name{Co5d5}
\alias{Co5d5}
\title{A data generating model with two clusters}
\usage{
Co5d5(n)
}
\arguments{
\item{n}{number of observations in the dataset. These are equally divided between clusters.}
}
\value{
Returns a list having two components described as follows;
 \describe{
 \item{n}{number of observations.}
  \item{K}{number of clusters.}
   \item{cluster size}{number of observations in each cluster.}
 \item{data}{data set generated from the model.}
\item{truelab}{clustering label vector correspounding to the known data generating model.}
}
}
\description{
Generates data set consists of five correlated clusters, one each from Gaussian distributions.
}
\details{
The data set has five correlated dimensions.
}
\examples{
dmat <- Co5d5(1000)
plot2d(dmat$data, dmat$truelab)
pairplots(dmat$data, dmat$truelab) 

}
\references{
Batool F. (2019). Optimum average silhouette width clustering.  
\emph{PhD Thesis}, University College London.
}
\author{
Fatima Batool \email{ucakfba@ucl.ac.uk}
}
