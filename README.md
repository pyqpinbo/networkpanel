# networkpanel
This package can identify the latent group structure for panel data models via community detection

# Example code
K  <- 3      # number of groups/classes
n  <- 60    # cross section dimension
T  <- 150    # time series dimension
bw <- 0.5*T^(-1/5)   # bandwidth for NW estimators of m_i
N  <- 100    # grid length for the estimation of the functions m_i, the functions are estimated on the grid 1/N, 2/N, ..., 1

# Simulation of data
# ----------------------------------------------
# generate regressors X
X.mat <- runif(n*T)
X.mat <- matrix(X.mat,ncol=n)
X <- list()
for(i in 1:n)
  X[[i]] <- X.mat[,i]
# generate error terms
sd.eps  <- sqrt(1)
eps.mat <- rnorm(n*T,mean=0,sd=sd.eps)
eps.mat <- matrix(eps.mat,ncol=n)
# generate g-functions
theta <- function(u){(1-u^2)^4 *((sign(1-u^2)+1)/2)}
g2 <- function(u){2.5*theta((u-0.75)/0.8) - 0.75}
g3 <- function(u){1.75*atan(5*(u-0.6)) + 0.75}
g1 <- function(u){sin(2*pi*u)}
# define group lengths
n2 <- ceiling(n/3)
n3 <- ceiling(n/3)
n1 <- n-n2-n3
nk.vec <- c(0,n1,n1+n2,n1+n2+n3)
# generate Y-observations
Y.mat <- matrix(0,ncol=n,nrow=T)
for(k in 1:K)
{  for(i in (nk.vec[k]+1):nk.vec[k+1])
{  g.mat <- cbind(g1(X.mat[,i]),g2(X.mat[,i]),g3(X.mat[,i]))
Y.mat[,i] <- g.mat[,k] + eps.mat[,i]
}
}

Y <- list()
for(i in 1:n)
  Y[[i]] <- Y.mat[,i]

# apply the proposed procedure to identify the latent group structure
library("networkpanel")
library("compiler")
library("parallel")
library("TSclust")
library("scales")

# Nadaraya-Watson estimates of the subject-specific functions m_i
m.hat <- matrix(0,ncol=n,nrow=N)   # Nadaraya-Watson estimates of the functions m_i; 
for(i in 1:n)
  m.hat[,i] <- m.i.hat(X[[i]],Y[[i]],bw,N)

ts_list = lapply(seq_len(ncol(m.hat)), function(i) m.hat[,i])
matplot(do.call(cbind, ts_list), t='l', lty=1, xlab="x", ylab="y",main="Initial estimator")

# Calculate the distance matrix and normalize it
distance_matrix = dist.parallel(tsList = ts_list, distFunc = tsdiss.euclidean, cores = 1)
distance_matrix = dist_normalize(distance_matrix)
heatmap(distance_matrix, main = "Distance matrix")

# Network construction 
net = network_build(distance_matrix = distance_matrix, percentile_value = 0.3)
net_layout = layout_components(net)

# plot the network
plot(net, layout=net_layout)

# community detection method to identify the latent group structure 
communities = cluster_louvain(net)
plot(communities, net, layout=net_layout)


