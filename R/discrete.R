
#' Get reduced density matrix
#'
#'@param bipartitite.qubits tensor of bipartite systems
#'@param keep.dim dimension to keep (default: 1)
#'@return reduced density matrix
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
#'@param red.dm reduced density matrix
#'@return Schmidt modes, including the eigenvalues, and eigenvectors of both subsystems of the modes
#'@export
schmidt.decompose<- function(bipartite.qubits) {
  mindim<- min(dim(qubits))
  red.dm<- reduced.denmat(bipartite.qubits)
  decomposed<- eigen(red.dm, TRUE)
  vecmat2<- bipartite.qubits %*% solve(decomposed$vectors)
  modes<- lapply(1:mindim,
                 function(i) list(eigenvalue=decomposed$values[i],
                                  sys1vector=decomposed$vectors[,i],
                                  sys2vector=vecmat2[,i]/sqrt(decomposed$values[i]))
                 )

  modes
}
