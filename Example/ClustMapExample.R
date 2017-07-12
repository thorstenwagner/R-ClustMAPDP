library(ClustMAPDP);
library(expm);
library(matrixStats);
set.seed(123)
#repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

D <- 2; # Data dimensionality
K <- 5; # Number of clusters

#Generate random cluster locations and covariance matrices
manyZeros <- rep(0,D*D*K);
Sigma <- array(manyZeros,c(D,D,K));
Mu <- matrix(0,D,K);

for(k in 1:K){
  s <- matrix(rnorm(D*D),D);
  s <- t(s)%*%s;
  Sigma[,,k]<-0.15*s;
  Mu[,k] <- 2.8*matrix(rnorm(D),D)
}
# Generate random cluster mixture weights
Pi <- matrix(runif(K),1);
Pi <- Pi/sum(Pi);

N <- 4000;

# Generate categorical data for cluster assignments z
Z <- sample(K,N,replace = TRUE,prob=Pi);

#Generate multivariate Gaussian data for X
X <- matrix(0,D,N);
for(k in 1:K){
  i <- (Z==k);
  M <- sum(i);
  X[,i] <- sqrtm(Sigma[,,k])%*%matrix(rnorm(D*M),D)+repmat(Mu[,k],1,M);
}
i <- sample(1:N,N,replace=FALSE);
Z <- Z[i];
X <- X[,i];
# Set up Normal-Wishart MAP-DP prior parameters
N0 <- 1;                           # Prior count (concentration parameter)
m0 <- rowMeans(X);                     # Normal-Wishart prior mean
a0 <- 10;                            # Normal-Wishart prior scale
c0 <- 10/N;                          # Normal-Wishart prior degrees of freedom
B0 <- diag(1./(0.05*rowVars(X)));    # Normal-Wishart prior precision
start.time <- Sys.time()
Rprof(tmp <- tempfile())
r <- clustMapDP(X,N0,m0,a0,c0,B0);
Rprof()
summaryRprof(tmp)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

par(mfrow=c(1,3));
plot(X[1,],X[2,],col=Z, main = "Input")
plot(X[1,],X[2,],col=r$z, main = paste("MAP-DP (",r$K,")"))
km <- kmeans(x=t(X),c=r$K,iter.max = 20)
plot(X[1,],X[2,],col=km$cluster, main = "kMeans")

