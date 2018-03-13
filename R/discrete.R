
#' Get reduced density matrix
#'
#'@param bipartitite.qubits tensor of bipartite systems
#'@param keep.dim dimension to keep (default: 1)
#'@export
reduced.denmat<- function(bipartite.qubits, keep.dim=1) {
  if (keep.dim==2) bipartite.qubits<- Conj(t(bipartite.qubits))
  keepdimnum<- dim(bipartite.qubits)[1]
  tracedimnum<- dim(bipartite.qubits)[2]
  rdenmat<- matrix(c(0)*(keepdimnum*keepdimnum), byrow=TRUE, nrow=keepdimnum)
  for (i in 1:keepdimnum) {
    for (j in 1:keepdimnum) {
      for (k in 1:tracedimnum) {
        rdenmat[i, j] = rdenmat[i, j] + bipartite.qubits[i, k]*Conj(bipartite.qubits[j, k])
      }
    }
  }
  if (keep.dim==2) Conj(t(rdenmat)) else rdenmat
}


