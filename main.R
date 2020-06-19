#' Distance matrix normalization
#'
#' This function normalize the distance matrix
#' @param dist The original distance matrix
#' @keywords distance matrix
#' @export
#' @examples
dist_normalize <- function(dist) {
    distNorm = matrix(0, nrow(dist), ncol(dist))
    d = dist[upper.tri(dist)]
    d = rescale(d)
    distNorm[upper.tri(distNorm)] = d
    distNorm = distNorm + t(distNorm)
    distNorm
}

#' Pearson correlation between time series
#'
#' This function calcualte pearson correlation between time series
#' @param ts1 The first time series
#' @param ts2 The second time series
#' @keywords pearson correlation
#' @export
#' @examples

tsdiss.correlation <- function(ts1, ts2) {
    CorDistance(ts1, ts2)
}

#' Euclidean L2 distance between time series
#'
#' This function calcualte L2 distance between time series
#' @param ts1 The first time series
#' @param ts2 The second time series
#' @keywords L2 distance
#' @export
#' @examples

tsdiss.euclidean <- function(ts1, ts2) {
    diss.EUCL(ts1, ts2)
}

#' Manhattan L1 distance between time series
#'
#' This function calcualte L1 distance between time series
#' @param ts1 The first time series
#' @param ts2 The second time series
#' @keywords L1 distance
#' @export
#' @examples

tsdiss.manhattan <- function(ts1, ts2) {
    sum(abs(ts1-ts2))
}

#' Infinite Norm between time series
#'
#' This function calcualte infinite norm between time series
#' @param ts1 The first time series
#' @param ts2 The second time series
#' @keywords L infinity distance
#' @export
#' @examples

tsdiss.infiniteNorm <- function(ts1, ts2) {
    max(abs(ts1-ts2))
}

#' Time series clustering using community detection
#'
#' This function conduct time series clustering using community detection
#' @param dist Distance matrix between every pair of time series.
#' @param epsilon Parameter for the network construction method.
#' @param communityDetectionFunc Community detection function. The igraph package has some community detection algorithms implemented.
#' @keywords community detection
#' @export
#' @examples

ts.community.detection <- function(dist,epsilon, communityDetectionFunc=cluster_louvain){
  net = net.epsilon.create(dist, epsilon)
  communities = communityDetectionFunc(net)
  communities=membership(communities)
  #return(communities)
}

#' Network construction using epsilon threshold
#'
#' This function creats network using threshold epsilon
#' @param dist Distance matrix between every pair of time series.
#' @param epsilon Parameter for the network construction method.
#' @keywords network construction
#' @export
#' @examples

net.epsilon.create <- function(dist, epsilon) {
  n = matrix(0, ncol(dist), nrow(dist))
  n[dist < epsilon] = 1;
  graph.adjacency(n, mode="undirected", diag=F)
}

#' Epanechnikov kernel function
#' @param x Variable
#' @keywords Epanechnikov kernel
#' @export
#' @examples

epan <- function(x)
{  return(0.75*(1-x^2)*((sign(1-x^2)+1)/2))}


#' NW estimator of nonparametric function
#' @param X indepent variable with dimension T*1
#' @param Y dependent variable with dimension T*1
#' @param bw bandwidth
#' @param N number of grid points
#' @keywords NW estimator
#' @export
#' @examples
m.i.hat <- function(X,Y,bw,N)
{  m.vec <- rep(0,N)
  for(j in 1:N)
  {  rh <- sum(epan((X-j/N)/bw) * Y)
     fh <- sum(epan((X-j/N)/bw))
     m.vec[j] <- rh/fh
  }
return(m.vec)
}


#' Parallel distance matrix
#' @param tsList indepent variable with dimension T*1
#' @param distFunc dependent variable with dimension T*1
#' @param cores bandwidt
#' @keywords Parallel distance matrix
#' @export
#' @examples

dist.parallel <- function(tsList, distFunc=tsdiss.euclidean, cores=2) {
  distFuncCompiled <- cmpfun(distFunc)
  tsListLength = length(tsList)
  combs = combn(tsListLength, 2)
  d = mcmapply(dist.parallel2.compute, x=combs[1,], y=combs[2,],
               MoreArgs=list(tsList=tsList, distFunc=distFuncCompiled), mc.cores=cores)
  dist = matrix(0, tsListLength, tsListLength)
  dist[lower.tri(dist)] = d
  dist = as.matrix(as.dist(dist))
  return(dist)
}

#' calculate distance matrix
#' @param x the first time series
#' @param y the second time series
#' @param tsList independent variable with dimension T*1
#' @param distFunc dependent variable with dimension T*1
#' @keywords Parallel distance matrix
#' @export
#' @examples
dist.parallel2.compute <- function(x, y, tsList, distFunc) {
  distFunc(tsList[[x]], tsList[[y]])
}

