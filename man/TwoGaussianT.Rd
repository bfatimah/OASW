% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoGaussianT.R
\name{TwoGaussianT}
\alias{TwoGaussianT}
\title{A data generating model with three clusters}
\usage{
TwoGaussianT(n)
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
Generates data set consists of three clusters, two clusters from Gaussian distributions 
which are compact and closely located to each other but far from the third cluster generated from the Student's t 
distribution with wider spread.
}
\details{
The data set has two iid dimensions.
}
\examples{
dmat <- TwoGaussianT(1000)
plot2d(dmat$data, dmat$truelab)
}
\references{
Batool F. (2019). Optimum average silhouette width clustering.  
\emph{PhD Thesis}, University College London.
}
\author{
Fatima Batool \email{ucakfba@ucl.ac.uk}
}
