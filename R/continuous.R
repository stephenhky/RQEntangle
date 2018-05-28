
#' Interpolate values of functions.
#'
#'@param xarr a vector of x (sorted)
#'@param yarr a vector of y
#'@param x given value of x
#'@return interpolated value of y
continuous.function.interpolate<- function(xarr, yarr, x) {
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


#' Lambda function of the interpolated continous function.
#'
#'@param xarr a vector of x (sorted)
#'@param yarr a vector of y
#'@return interpolated lambda function
interpolated.continuous.function<- function(xarr, yarr) {
  function(x) mapply(function(xi) continuous.function.interpolate(xarr, yarr, xi), x)
}


#' Making a discretized tensor for a continuous function
#'
#'@param bifunc bipartitite continuous wavefunction
#'@param x1lo lower limit of \code{x1}
#'@param x1hi upper limit of \code{x1}
#'@param x2lo lower limit of \code{x2}
#'@param x2hi upper limit of \code{x2}
#'@param nbx1 number of discretized x1 (default: 100)
#'@param nbx2 number of discretized x2 (default: 100)
#'@return discretized tensor for Schmidt decomposition
#'@importFrom itertools ihasNext hasNext product
#'@importFrom iterators nextElem
discretize.continuous.bipartitefunc<- function(bifunc, x1lo, x1hi, x2lo, x2hi, nbx1=100, nbx2=100) {
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

#' Perform a continuous Schmidt decomposition
#'
#'@param bifunc bipartitite continuous wavefunction
#'@param x1lo lower limit of \code{x1}
#'@param x1hi upper limit of \code{x1}
#'@param x2lo lower limit of \code{x2}
#'@param x2hi upper limit of \code{x2}
#'@param nbx1 number of discretized x1 (default: 100)
#'@param nbx2 number of discretized x2 (default: 100)
#'@param keep number of Schmidt modes to keep (default: minimum of 10, \code{nbx1}, and \code{nbx2})
#'@return Schmidt modes, including the eigenvalues, and the lambda interpolated function of the Schmidt modes
#'
#'@examples
#'coupled.harm.fcn<- function(x1,x2) exp(-((0.5*(x1+x2))**2))*exp(-(x1-x2)**2)*sqrt(2./pi)
#'continuous.schmidt.decompose(coupled.harm.fcn, -10, 10, -10, 10)
#'
#'@export
continuous.schmidt.decompose<- function(bifunc, x1lo, x1hi, x2lo, x2hi, nbx1=100, nbx2=100, keep=min(10, nbx1, nbx2)) {
  tensor<- discretize.continuous.bipartitefunc(bifunc, x1lo, x1hi, x2lo, x2hi, nbx1=nbx1, nbx2=nbx2)
  dis.decomposition<- schmidt.decompose(tensor)

  dx1<- (x1hi-x1lo)/(nbx1-1)
  dx2<- (x2hi-x2lo)/(nbx2-1)
  x1array<- (1:nbx1)*dx1+x1lo
  x2array<- (1:nbx2)*dx2+x2lo

  sum.eigvals<- sum(mapply(function(decomposition) decomposition$eigenvalue, dis.decomposition))
  lapply(dis.decomposition,
         function(decomposition) {
           sqnorm1<- sum(decomposition$sys1vector^2) * dx1
           sqnorm2<- sum(decomposition$sys2vector^2) * dx2
           list(eigenvalue=decomposition$eigenvalue / sum.eigvals,
                sys1eigfcn=interpolated.continuous.function(x1array, decomposition$sys1vector / sqrt(sqnorm1)),
                sys2eigfcn=interpolated.continuous.function(x2array, decomposition$sys2vector / sqrt(sqnorm2))
           )
         })
}

