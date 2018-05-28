
#' Calculate the entanglement entropy given the calculate Schmidt modes.
#'
#'@param modes Schmidt modes
#'@return entanglement entropy
#'
#'@examples
#'singlet<- matrix(c(0, sqrt(0.7), sqrt(0.3), 0), byrow = TRUE, nrow = 2)
#'modes<- schmidt.decompose(singlet)
#'entanglement.entropy(modes)
#'
#'@export
entanglement.entropy<- function(modes) {
  eigenvalues<- mapply(function(mode) mode$eigenvalue, modes)
  eigenvalues<- eigenvalues[ eigenvalues>0]
  -sum(eigenvalues*log(eigenvalues))
}


#' Calculate the participation ratio given the calculate Schmidt modes.
#'
#'@param modes Schmidt modes
#'@return participation ratio
#'
#'@examples
#'singlet<- matrix(c(0, sqrt(0.7), sqrt(0.3), 0), byrow = TRUE, nrow = 2)
#'modes<- schmidt.decompose(singlet)
#'participation.ratio(modes)
#'
#'@export
participation.ratio<- function(modes) {
  eigenvalues<- mapply(function(mode) mode$eigenvalue, modes)
  eigenvalues<- eigenvalues[ eigenvalues>0]
  1/sum(eigenvalues*eigenvalues)
}


#' Calculate the negativity given the calculate Schmidt modes.
#'
#'@param modes Schmidt modes
#'@return negativity
#'
#'@examples
#'singlet<- matrix(c(0, sqrt(0.7), sqrt(0.3), 0), byrow = TRUE, nrow = 2)
#'modes<- schmidt.decompose(singlet)
#'negativity(modes)
#'
#'@export
negativity<- function(modes) {
  eigenvalues<- mapply(function(mode) mode$eigenvalue, modes)
  0.5*(sum(abs(eigenvalues))-1)
}
