
#' Calculate the entanglement entropy given the calculate Schmidt modes.
#'
#'@modes calculated Schmidt modes
#'@return entanglement entropy
#'@export
entanglement.entropy<- function(modes) {
  eigenvalues<- mapply(function(mode) mode$eigenvalue, modes)
  eigenvalues<- eigenvalues[ eigenvalues>0]
  -sum(eigenvalues*log(eigenvalues))
}


#' Calculate the participation ratio given the calculate Schmidt modes.
#'
#'@modes calculated Schmidt modes
#'@return participation ratio
#'@export
participation.ratio<- function(modes) {
  eigenvalues<- mapply(function(mode) mode$eigenvalue, modes)
  eigenvalues<- eigenvalues[ eigenvalues>0]
  1/sum(eigenvalues*eigenvalues)
}


#' Calculate the negativity given the calculate Schmidt modes.
#'
#'@modes calculated Schmidt modes
#'@return negativity
#'@export
negativity<- function(modes) {
  eigenvalues<- mapply(function(mode) mode$eigenvalue, modes)
  0.5*(sum(abs(eigenvalues))-1)
}
