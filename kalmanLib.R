#
# This is a set of functions which implement many of the necessary calculations
# for a Kalman filter with arbitrary state.  Matrix operators and dimensions of 
# input matricies are used.
#
#

#
# Identity matrix
#
eye <- function(n=3){ diag(x=1,ncol=n,nrow=n)}

#
# Create a column vector of zeros
#
colVector <- function(x=0,n=3){ matrix(x,ncol=1,nrow=n)}
#
# Create a square matrix of zeros
#
sqMatrix <- function(x=0,n=3){ matrix(x,ncol=n,nrow=n)}
#
# This function returns a time-translation operator that 
# implements newtonian physics
#
newtonTimeTranslation <- function(n=6) {
  matrix(c(c(1,0,0,0,0,0),
           c(0,1,0,0,0,0),
           c(1,0,1,0,0,0),
           c(0,1,0,1,0,0),
           c(0.5,0,1,0,1,0),
           c(0,0.5,0,1,0,1)),nrow=n,ncol=n) # time-evolution operator
}
  


#
# Perform prediction step -- like a time-translation operator in QM
# (use physics or other rules about how the system propagates forward
# in time).  
# 
# Fm is the time-translation operator matrix.  
# G is the "control-input" model, which is just a way of tranforming
# a known input (u) into state-space (x).
#
pred <- function(x){ Fm%*%x + G%*%u}
#
# This function takes its name from the intended use (to measure the 
# state of the system after the prediction step), but what it is really
# doing is applying a matrix (H) to transform from state-space (x) to 
# observation-space (z).  Note that H is used to calculate the "Kalman gain", 
# which is more than just a gain term (though it is that also).  The Kalman 
# gain transforms residuals in observation-space (z) to state-space (x).
#
measurePred <- function(x){ H%*%x }
#
# Perform the residual calculation step.  This function probably does too
# much at the moment: Currently, the function calculates the residual
# between the simulated actual observation-state of the system (which the Kalman filter
# toy pretends to not know) and the observed-space state of the system (z).
# 
# Random noise is also added to the observation in this function (note this
# is noise in the observation itself)
#
# Returned value is the measurement residual (v) in observation-space.
residual <- function(z,k){ 
  res <- H%*%physMod(k) - z 
  noise <- H %*% matrix(rnorm(n=3,sd=0.1),nrow=dim(z)[1],ncol=dim(z)[2])
  return(res + noise)
}

#
# This function performs the update step of the Kalman filter by mixing
# the predicted state (x) with the measurement residual (v).  Note that 
# x is in state-space and v is in observation-space, so the Kalman gain
# trasformation matrix (W) does double duty: It transforms v into state-space
# and also acts as a gain term.  
# 
# Note that this tranformation is what allows measurements of only (e.g.)
# position to affect the system's acceleration state.
#
# At the high-gain (perfect measurements) 
# limit, the predicted state is reset to the observed state.  At the
# low gain limit (very noisy measurements), the updated predicted state
# is almost unaffected by the measurement.
#
update <- function(x,v,W){ x + W%*%v}

#
# The gain function uses the state prediction covariance and measurement
# prediction covariance to compute the Kalman gain.
#
# Fm is the time-translation operator
# P is the a-priori state prediction/estimate covariance
# Pko is the post-posteriori state prediction/estimate covariance
# Q is the covariance of the "process noise" (noise which internally
# perturbs the state of the system)
# H is the state-space ==> observation-space tranformation operator
# R is the covariance of the "observation noise"
# Sko is the post-posteriori measurement covariance
#
# More info:
#http://biorobotics.ri.cmu.edu/papers/sbp_papers/integrated3/kleeman_kalman_basics.pdf
#
gain <- function(P){
  Pko <- Fm%*%P%*%t(Fm) + Q
  Sko <- H%*%Pko%*%t(H) + R
  return(Pko%*%t(H)%*%solve(Sko))
}

#
# This function updates the estimate of the state prediction/estimate
# covariance from the filter gain (W), time-translation operator (Fm),
# process noise (Q), observation noise (R) and observation transformation (H)
#
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
