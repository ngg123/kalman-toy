eye <- function(n=3){ diag(x=1,ncol=n,nrow=n)}

colVector <- function(x=0,n=3){ matrix(x,ncol=1,nrow=n)}
sqMatrix <- function(x=0,n=3){ matrix(x,ncol=n,nrow=n)}

pred <- function(x){ Fm%*%x + G%*%u}

measurePred <- function(x){ H%*%x }

residual <- function(z,k){ H%*%physMod(k)-z + H%*%matrix(rnorm(n=3,sd=0.1),nrow=dim(z)[1],ncol=dim(z)[2])}

update <- function(x,v,W){ x + W%*%v}

gain <- function(P){
  Pko <- Fm%*%P%*%t(Fm) + Q
  Sko <- H%*%Pko%*%t(H) + R
  return(Pko%*%t(H)%*%solve(Sko))
}

covP <- function(W){
  #   Pko <- Fm%*%P%*%t(Fm) + Q
  #   Sko <- H%*%Pko%*%t(H) + R
  #   return(Pko - W%*%Sko%*%t(W))
  
  # "Joseph's form" covariance update to 
  # avoid catastrophic cancellation
  Pko <- Fm%*%P%*%t(Fm) + Q
  SUB <- diag(dim(W)[1]) - W%*%H
  return(SUB%*%Pko%*%t(SUB) + W%*%R%*%t(W))
}
