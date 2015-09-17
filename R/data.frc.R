data.frc <- function(data.in,method=c("crost","crost.ma","tsb","sexsm","imapa"),...){
# Wrapper to forecasts data.frames with a single call
# 
# Inputs:
#   data.in     Data frame with time series. This can also be a matrix or array with each 
#               column being a different time series. 
#   method      Which method to use for forecasting: 
#                 "crost", "crost.ma", "tsb", "sexsm", "imapa"
#   ...         Additional inputs to pass to forecasting functions. See individual function 
#               documentation for options. 
#
# Outputs:
#   frc.out     Data frame containing forecasts for all time series.
#   out         List with detailed output per series. To access individual outputs of the list
#               use: sapply(out, get, x="element"), where "element" could be for example "frc.in". 
#
# Nikolaos Kourentzes, 2015 <nikolaos@kourentzes.com>

  # Get defaults
  method <- method[1]
  
  # Get number of columns and class of data.in
  # k <- ncol(data.in)
  # data.cl <- class(data.in)
  data.names <- colnames(data.in)
  
  # Check is selected model is implemented
  allow.method <- c("crost","crost.ma","tsb","sexsm","imapa")
  if (!(method %in% allow.method)){
    stop(paste(method,"not a permitted forecasting method."))
  }
  
  # Fit models and produce forecasts
  switch(method,
         "crost" = {out <- apply(data.in,2,crost,...)},
         "crost.ma" = {out <- apply(data.in,2,crost.ma,...)},
         "tsb" = {out <- apply(data.in,2,tsb,...)},
         "sexsm" = {out <- apply(data.in,2,sexsm,...)},
         "imapa" = {out <- apply(data.in,2,imapa,...)}
         )
  
  # Now get from res all forecasts and return them to a data.frame
  frc.out <- sapply(out, get, x="frc.out")
  colnames(frc.out) <- data.names
  frc.out <- data.frame(frc.out)
  
  return(list(frc.out = frc.out, out = out))

}