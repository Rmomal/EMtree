##################################################################
# Miscelaneous functions
##################################################################



#' Makes a symmetric matrix from the vector made of its upper triangular part
#'
#' @param A.vec Vector obtain with ToVec or the upper.tri() function.
#'
#' @return The symmetric matrix with A.vec as off-diagonal terms.
#' @export
ToSym <- function(A.vec){
   n = (1+sqrt(1+8*length(A.vec)))/2
   A.mat = matrix(0, n, n)
   A.mat[upper.tri(A.mat)] = A.vec
   A.mat = A.mat + t(A.mat)
   return(A.mat)
}


#' Makes a vector from the upper triangular part of a symmetric matrix
#'
#' @param A.mat A symmetric matrix.
#'
#' @return The vector from the upper triangular part of A.mat.
#' @export
ToVec <- function(A.mat){
   return(A.mat[upper.tri(A.mat)])
}


