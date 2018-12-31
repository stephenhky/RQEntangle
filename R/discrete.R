
#' Get reduced density matrix
#'
#'@param bipartite.qubits tensor of bipartite systems
#'@param keep.dim dimension to keep (default: 1)
#'@return reduced density matrix
#'
#'@examples
#'singlet<- matrix(c(0, sqrt(0.7), sqrt(0.3), 0), byrow = TRUE, nrow = 2)
#'reduced.denmat(singlet)
#'
#'@export
reduced.denmat<- function(bipartite.qubits, keep.dim=1) {
  if (keep.dim==2) bipartite.qubits<- Conj(t(bipartite.qubits))
  keepdimnum<- dim(bipartite.qubits)[1]
  tracedimnum<- dim(bipartite.qubits)[2]
  rdenmat<- matrix(rep(0, keepdimnum*keepdimnum), byrow=TRUE, nrow=keepdimnum)
  for (i in 1:keepdimnum) {
    for (j in 1:keepdimnum) {
      for (k in 1:tracedimnum) {
        rdenmat[i, j] = rdenmat[i, j] + bipartite.qubits[i, k]*Conj(bipartite.qubits[j, k])
      }
    }
  }
  if (keep.dim==2) Conj(t(rdenmat)) else rdenmat
}


#' Perform Schmidt decomposition
#'
#'@param bipartite.qubits tensor of bipartite systems
#'@return Schmidt modes, including the eigenvalues, and eigenvectors of both subsystems of the modes
#'
#'@examples
#'singlet<- matrix(c(0, sqrt(0.7), sqrt(0.3), 0), byrow = TRUE, nrow = 2)
#'schmidt.decompose(singlet)
#'
#'@export
schmidt.decompose<- function(bipartite.qubits) {
  svdres<- svd(bipartite.qubits)
  nb_modes<- length(svdres$d)
  modes<- lapply(1:nb_modes,
                 function(i) list(eigenvalue=svdres$d[i]*svdres$d[i],
                                  sys1vector=svdres$u[,i],
                                  sys2vector=svdres$v[,i])
                 )

  lapply(order(mapply(function(mode) mode$eigenvalue, modes), decreasing = TRUE), function(i) modes[[i]])
}
