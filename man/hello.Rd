\name{Entropy Weighted Dirichlet Process Means}
\alias{EWDP-Means}
\title{EWDPmeans}
\usage{
hello()
}
\description{
Prints 'This package is a clustering technique that can provide partitioning of data in an unsupervised way. The additional advantage of using this method is that the number of clusters are not required to known a-priori. This method is useful in clustering high-dimensional data since the entropy regularizer on feature weights provides a method of feature selection.
The input requires a properly centred and scaled data-matrix in data.matrix form and additional parameters lambda_w, lambda_k and tmax. The output provides the number of clusters, feature weights and the cluster labels for each data-point.'
}
\examples{
  #Loading the data
  data(iris)
  X=iris
  X=data.matrix(X)
  toss=X[,5]
  X=X[,-5]
  #Centre and scale
  p=dim(X)[2]
  for(l in 1:p){
    X[,l]=(X[,l]-mean(X[,l]))/sd(X[,l])
  }
  #Run the code
  l=EWDPmeans(X,1,1,30)
  #Number of clusters
  l[[1]]
  #Feature weights
  l[[2]]
  #Cluster labels
  l[[3]]
}
