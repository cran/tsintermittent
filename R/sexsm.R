sexsm <- function(data,h=10,w=NULL,init=c("mean","naive"),cost=c("mar","msr","mae","mse"),
                  init.opt=c(TRUE,FALSE),outplot=c(FALSE,TRUE),opt.on=c(FALSE,TRUE)){
# Simple Exponential Smoothing
#
# Inputs:
#   data        Intermittent demand time series.
#   h           Forecast horizon.
#   w           Smoothing parameter. If w == NULL then parameter is optimised.
#   init        Initial values for demand and intervals. This can be:
#                 x       - Numeric value for the initial level;
#                 "naive" - Initial value is a naive forecast;
#                 "mean"  - Initial value is equal to the average of data.
#   cost        Cost function used for optimisation
#                 "mar" - Mean absolute rate
#                 "msr" - Mean squared rate
#                 "mae" - Mean absolute error
#                 "mse" - Mean squared error
#   init.opt    If parameters are optimised and init.opt==TRUE then initial values are 
#               optimised as well. 
#   outplot     If TRUE a plot of the forecast is provided.
#   opt.on      This is meant to use only by the optimisation function. When opt.on is 
#               TRUE then no checks on inputs are performed. 
#
# Outputs:
#   model       Type of model fitted.
#   frc.in      In-sample demand. 
#   frc.out     Out-of-sample demand.
#   alpha       Smoothing parameter.
#   initial     Initialisation value.
#
# Example:
#   sexsm(ts.data1,outplot=TRUE)
#
# Notes:
# Optimisation of the method described in:
# N. Kourentzes, 2014, On intermittent demand model optimisation and selection, 
# International Journal of Production Economics, 156: 180-190. 
# http://dx.doi.org/10.1016/j.ijpe.2014.06.007
# http://kourentzes.com/forecasting/2014/06/11/on-intermittent-demand-model-optimisation-and-selection/
#
# Nikolaos Kourentzes, 2014 <nikolaos@kourentzes.com>
  
  # Defaults
  w <- w[1]
  cost <- tolower(cost[1])
  init.opt <- init.opt[1]
  outplot <- outplot[1]
  opt.on <- opt.on[1]
  if (!is.numeric(init)){
    init <- init[1]
  } 
  
  n <- length(data)
  
  # Initialise
  if (!is.numeric(init)){
    if (init=="mean"){
      init <- mean(data)
    } else {
      init <- data[1]
    }
  }
  
  # Optimise parameters if requested
  if (is.null(w)){
    wopt <- sexsm.opt(data,cost,init,init.opt)
    w <- wopt$w
    init <- wopt$init
  }
  
  # Pre-allocate memory
  fit <- vector("numeric",n+1)
  
  # Assign initial values and parameters
  if (opt.on == FALSE){
    fit[1] <- init[1]
  }
  
  # Fit model
  for (i in 2:(n+1)){
    fit[i] <- w*data[i-1] + (1-w)*fit[i-1]
  }
  
  # Calculate in-sample demand rate
  frc.in <- fit[1:n]
  
  # Forecast out-of-sample demand rate
  if (h>0){
    frc.out <- rep(fit[n+1],h)
  } else {
    frc.out = NULL
  }
  
  # Plot
  if (outplot==TRUE){
    plot(1:n,data,type="l",xlim=c(1,(n+h)),xlab="Period",ylab="",
         xaxs="i",yaxs="i",ylim=c(0,max(data)*1.1))
    lines(which(data>0),data[data>0],type="p",pch=20)
    lines(1:n,frc.in,col="red")
    lines((n+1):(n+h),frc.out,col="red",lwd=2)
  }
  
  return(list(model="sexsm",frc.in=frc.in,frc.out=frc.out,
              alpha=w,initial=init))

}

#-------------------------------------------------
sexsm.opt <- function(data,cost=c("mar","msr","mae","mse"),
                      init,init.opt=c(TRUE,FALSE)){
# Optimisation function for SES
  
  cost <- cost[1]
  init.opt <- init.opt[1]  
  
  if (init.opt == TRUE){
    w <- c(0.5,init[1])
    lbound <- c(0,0)
    ubound <- c(1,max(data))
    wopt <- optim(par=w,sexsm.cost,method="Nelder-Mead",data=data,cost=cost,
                  init=init,init.opt=init.opt,lbound=lbound,ubound=ubound,
                  control=list(maxit = 2000))$par
  } else {
    w <- 0.5
    lbound <- 0
    ubound <- 1
    wopt <- optim(par=w,sexsm.cost,method="Brent",data=data,cost=cost,
                  init=init,init.opt=init.opt,lower=lbound,upper=ubound,
                  control=list(maxit = 2000))$par    
    wopt <- c(wopt,init)
  }
  
  return(list(w=wopt[1],init=wopt[2]))
  
}

#-------------------------------------------------
sexsm.cost <- function(w,data,cost,init,init.opt,lbound,ubound){
# Cost functions for SES
  
  if (init.opt==TRUE){
    frc.in <- sexsm(data=data,w=w[1],h=0,init=w[2],opt.on=TRUE)$frc.in
  } else {
    frc.in <- sexsm(data=data,w=w[1],h=0,init=init,opt.on=TRUE)$frc.in
  }
  
  if (cost == "mse"){
    E <- data - frc.in  
    E <- E[!is.na(E)]
    E <- mean(E^2)
  } else if(cost == "mae"){
    E <- data - frc.in  
    E <- E[!is.na(E)]
    E <- mean(abs(E))
  } else if(cost == "mar"){
    n <- length(data)
    temp <- cumsum(data)/(1:n)
    n <- ceiling(0.3*n)
    temp[1:n] <- temp[n]
    E <- abs(frc.in - temp)
    E <- E[!is.na(E)]
    E <- sum(E)
  } else if(cost == "msr"){
    n <- length(data)
    temp <- cumsum(data)/(1:n)
    n <- ceiling(0.3*n)
    temp[1:n] <- temp[n]
    E <- (frc.in - temp)^2
    E <- E[!is.na(E)]
    E <- sum(E)    
  }
  
  # Constrains
  if (init.opt==TRUE){
    for (i in 1:(1+1*init.opt)){
      if (!(w[i]>=lbound[i]) | !(w[i]<=ubound[i])){
        E <- 9*10^99
      }
    }
  }
  
  return(E)
  
}