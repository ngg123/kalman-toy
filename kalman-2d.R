

n <- 2*3 # two dimensions
wvar <- 2
vvar <- 2

#
# Define motion of object at k-th time step
#
#
# This is the model that the residual() function will use to calculate
# observations of the system. This is an arbitrary function that the Kalman
# filter is trying to follow.
#
xo <- matrix(c(0,0,80,0,-10,0),nrow=n,ncol=1) # initial state for phys engine
physMod <- function(k){ matrix(c(c(1,0,0,0,0,0),
                                 c(0,1,0,0,0,0),
                                 c(k,0,1,0,0,0),
                                 c(0,k,0,1,0,0),
                                 c(0.5*k^2,0,k,0,1,0),
                                 c(0,0.5*k^2,0,k,0,1)),
                               nrow=n,ncol=n)%*%c(xo[1:4],-9.8,0)}

source('./kalmanLib.R') 

x <- colVector(0,n) # system initial state
#
# Physics 
#
Fm <- newtonTimeTranslation()
G <- matrix(c(0.5,0,1,0,0,0),nrow=n,ncol=1) # control input -> state t-form
u <- 0 # control input

I <- eye(n) # identity matrix
W <- eye(n) # Kalman Gain (and obs->state t-form)
# state->obs t-form (H) turns x- and y-position state into observations
H <- matrix(c(c(1,0,0,0,0,0),
              c(0,1,0,0,0,0),
              c(0,0,0,0,0,0),
              c(0,0,0,0,0,0),
              c(0,0,0,0,0,0),
              c(0,0,0,0,0,0)
              ),nrow=n,ncol=n) 

#
# These noise covariance matricies are somewhat arbitrary for now
#
P <- eye(n) * 10 # state covariance "intuition"
#Q <- sqMatrix(x=0,n=n) # process noise covariance "all zeros"
Q <- diag(vvar,nrow=n,ncol=n)
R <- diag(wvar,nrow=n,ncol=n) # measurement noise covariance


runKal <- function(maxit=5){
  retval <- c()
  retx <- c()
  for(i in 1:maxit){
    x <- pred(x)
    retx <- rbind(retx,data.frame(i,cbind(t(x),t(physMod(i)))))
    zhat <- measurePred(x)
    v <- residual(zhat,i)
    #     print(v)
    retval <- rbind(retval,data.frame(i,t(physMod(i)-x)))
    x <- update(x,v,W)
    
    W <- gain(P)
    P <- covP(W)
  }
  colnam <- c('i','x','y','xdot','ydot','xddot','yddot',
              'x.p','y.p','xdot.p','ydot.p','xddot.p','yddot.p')
  names(retval) <- colnam[1:7]
  names(retx) <- colnam
    return(retx)
}


testKal <- function(){
  df <- runKal(maxit=50)
  plot(df$x)
  plot(df$y)
}