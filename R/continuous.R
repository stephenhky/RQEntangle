
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


#' Making a discretized tensor for a continuous function
#'
#'@param bifunc bipartitite continuous function
#'@param x1lo lower limit of \code{x1}
#'@param x1hi upper limit of \code{x1}
#'@param x2lo lower limit of \code{x2}
#'@param x2hi upper limit of \code{x2}
#'@param nbx1 number of discretized x1 (default: 100)
#'@param nbx2 number of discretized x2 (default: 100)
#'@export
discretize_continuous_bipartitefunc<- function(bifunc, x1lo, x1hi, x2lo, x2hi, nbx1=100, nbx2=100) {
  dx1<- (x1hi-x1lo)/(nbx1-1)
  dx2<- (x2hi-x2lo)/(nbx2-1)
  tensor<- matrix(rep(0, nbx1*nbx2), byrow = TRUE, nrow = nbx1)
  iterator<- ihasNext(product(i=1:nbx1, j=1:nbx2))
  while (hasNext(iterator)) {
    counter<- nextElem(iterator)
    tensor[counter$i, counter$j]<- bifunc(x1lo+(counter$i-1)*dx1, x2lo+(counter$j-1)*dx2)
  }
  tensor
}



