

#
# Define motion of object at k-th time step
#
xo <- matrix(c(0,100,-10),nrow=3,ncol=1) # initial state for phys engine
physMod <- function(k){ matrix(c(c(1,0,0),c(k,1,0),c(k^2/2,k,1)),nrow=3,ncol=3)%*%c(xo[1:2],-9.8)}

source('./kalmanLib.R')

x <- colVector(0,3) # state
Fm <- matrix(c(c(1,0,0),c(1,1,0),c(0.5,1,1)),nrow=3,ncol=3) # time-evolution operator
G <- matrix(c(0.5,1,0),nrow=3,ncol=1) # control input -> state t-form
u <- 0 # control input

I <- eye(3) # identity matrix
W <- eye(3) # Kalman Gain (and obs->state t-form)
H <- matrix(c(c(1,0,0),c(0,0,0),c(0,0,0)),nrow=3,ncol=3) # state->obs t-form
P <- eye(3) * 10 # state covariance "intuition
Q <- matrix(0,nrow=3,ncol=3) # process noise covariance "all zeros"
R <- diag(1,nrow=3,ncol=3) # measurement noise covariance


runKal <- function(maxit=5){
  retval <- c()
  for(i in 1:maxit){
    x <- pred(x)
    zhat <- measureP(x)
    v <- residual(zhat,i)
    #     print(v)
    retval <- rbind(retval,data.frame(i,t(physMod(i)-x)))
    x <- update(x,v,W)
    W <- gain(P)
    P <- covP(W)
  }
  names(retval) <- c('i','x','xdot','xddot')
  return(retval)
}