
#' Interpolate values of functions.
#'
#'@param xarr a vector of x (sorted)
#'@param yarr a vector of y
#'@param x given value of x
#'@return interpolated value of y
#'@export
continuous_function_interpolate<- function(xarr, yarr, x) {
  if (x==max(xarr)) {
    yarr[ length(yarr)]
  } else {
    # binaray search
    left<- 1
    right<- length(xarr)
    id<- floor((left+right)/2)
    while((id!=1) & (id!=length(xarr)) & !((x>=xarr[id]) & (x<xarr[id+1]))) {
      if (x>=xarr[id+1]) {
        left<- id+1
      } else if (x<xarr[id]) {
        right<- id-1
      }
      id<- floor((left+right)/2)
    }

    yarr[id]+(yarr[id+1]-yarr[id])/(xarr[id+1]-xarr[id])*(x-xarr[id])
  }
}


#' Lambda function of the interpolated continous function
#'
#'@param xarr a vector of x (sorted)
#'@param yarr a vector of y
#'@return interpolated lambda function
#'@export
interpolated_continuous_function<- function(xarr, yarr) {
  function(x) continuous_function_interpolate(xarr, yarr, x)
}
