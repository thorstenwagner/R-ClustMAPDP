library(matrixStats);
#' MAP-DP Clustering
#'
#'It implements the algorithm descriped by
#'Raykov YP, Boukouvalas A, Baig F, Little MA (2016)
#'What to Do When K-Means Clustering Fails: A Simple yet Principled Alternative Algorithm.
#'PLOS ONE 11(9): e0162259. https://doi.org/10.1371/journal.pone.0162259
#'
#' @param X DxN matrix of data
#' @param N0 prior count (DP concentration paramter)
#' @param m0 cluster prior mean
#' @param a0 cluster prior scale
#' @param c0 cluster prior degrees of freedom
#' @param B0 cluster prior precision (inverse covariance)
#' @export
clustMapDP <- function(X,N0,m0,a0,c0,B0) {
 D <- dim(X)[1];
 N <- dim(X)[2];
 epsilon <- 1e-6;

 #Initialization
 K <- 1;
 z <- matrix(1,N,1);
 Enew <- Inf;
 dE <- Inf;
 iter <- 0;
 E <- c();

 while (abs(dE) > epsilon) {

   Eold <- Enew;
   dik <- matrix(Inf,N,1);

   for(i in 1:N){

     dk <- matrix(Inf,K+1,1);

     Nki <- matrix(Inf,K+1,1);

     for(k in 1:K){

       zki <- (z==k);
       zki[i] <- FALSE;
       Nki[k] <- sum(zki);

       # Updates meaningless for Nki=0
       if(Nki[k]==0){
         next;
       }

       # Update NW cluster hyper parameters
        nwupdList <- nwupd(Nki[k],X[,zki],m0,a0,c0,B0);
        mki <- nwupdList$mki;
        aki <- nwupdList$aki;
        cki <- nwupdList$cki;
        Bki <- nwupdList$Bki;

       # Computer Student-t NLL, existing cluster
        help <-stnll(X[,i],mki,aki,cki,Bki,D);
        dk[k] <- help;
     }

     if(iter==0){
       Nki[1]<-1;
     }

     Nki[K+1] = N0;

     # Computer Student-t NLL, new cluster

     dk[K+1] <- stnll(X[,i],m0,a0,c0,B0,D);

     # Computer MAP Assignment
     v = dk-log(Nki);
     z[i] <- which.min(dk-log(Nki));
     dik[i] <- v[z[i]];

     #print(z[i])
     # Create new cluster if required

     if(z[i] == (K+1)){
       K <- K+1;
     }
   }

   # Remove any empty clusters and re-assign
   Knz = 0;
   for(k in 1:K){
     i <- (z==k);
     Nk <- sum(i);
     if (Nk>0){
       Knz <- Knz + 1;
       z[i] <- Knz;
     }
   }
   K <- Knz;
   Nk <- hist(z,1:K,plot=FALSE)$counts

   # Compute updated NLL
   Enew <- sum(dik)-K*log(N0)-sum(lgamma(Nk));
   dE <- Eold - Enew;
   iter <- iter +1;
   E <- c(E,Enew);
 }

 # Compute cluster centroids
 mu <- matrix(NaN,D,K);
 for(k in 1:K){
   xk <- X[,z==k];

   mu[,k] <- rowMeans(xk);
 }

 resultList <- list("mu"=mu,"z"=z,"K"=K,"E"=E);
 return(resultList);
}

#' Compute Stundent-t log likeligood
stnll <- function(x,m,a,c,B,D) {
  mu <- m;
  nu <- a-D+1;
  Lambda <- c*nu/(c+1)*B;
  x <- as.matrix(x);#matrix(x,nrow = 2,ncol = 1)
  nl <- (nu+D)/2*log(1+t(x-mu)%*%Lambda%*%(x-mu)/nu)-0.5*log(det(Lambda))+lgamma(nu/2)-lgamma((nu+D)/2)+D/2*log(nu*pi);

  return(as.numeric(nl))
}

repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

#' Update Normal-Wishart hyper parameters
nwupd <- function(Nki,xki,m0,a0,c0,B0) {
  xmki <- 0;
  if(!is.matrix(xki)){
    xmki <- mean(xki);
  }else{
    xmki <- rowMeans(xki);
  }
  if(is.matrix(xki)==FALSE){
    xki <- as.matrix(xki);
  }
  xmcki <- sweep(xki,1,repmat(xmki,1,Nki),"-"); #xki-repmat(xmki,1,Nki);

  Ski <- xmcki%*%t(xmcki);
  cki <- c0+Nki;
  mki <- (c0*m0+Nki*xmki)/cki;
  xm0cki <- xmki-m0;
  Bki <- solve(solve(B0)+Ski+c0*Nki/cki*xm0cki%*%t(xm0cki));
  aki <- a0+Nki;

  resultList <- list("mki"=mki,"aki"=aki,"cki"=cki,"Bki"=Bki);
  return(resultList);
}
