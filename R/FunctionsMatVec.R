##################################################################
# Miscelaneous functions
##################################################################
F_LowerTri2Sym <- function(A.mat, diag=F){
   # Makes a symmetric matrix from its lower tirangular
   A.sym = A.mat*lower.tri(A.mat) + t(A.mat*lower.tri(A.mat))
   if (diag){A.sym = A.sym + diag(diag(A.mat))}
   return(A.sym)
}

##################################################################
F_Vec2Sym <- function(A.vec){
   # Makes a symmetric matrix from the vector made of its lower tirangular part
   n = (1+sqrt(1+8*length(A.vec)))/2
   A.mat = matrix(0, n, n)
   A.mat[lower.tri(A.mat)] = A.vec
   A.mat = A.mat + t(A.mat)
   return(A.mat)
}

##################################################################
F_Sym2Vec <- function(A.mat){
   # Makes a vector from the lower triangular par of a symmetric matrix
   return(A.mat[lower.tri(A.mat)])
}
##################################################################
F_Array2Mat <- function(A.array){
  # Makes an n x pq matrix from a n x p x q array
  q = dim(A.array)[3]
  if (q==1){
    A.mat = A.array
  }else{
    A.mat = A.array[, , 1]
    invisible(sapply(2:q, function(k){A.mat <<- cbind(A.mat, A.array[, , k])}))
    }
  return(A.mat)
}
##################################################################
F_SymArray2Mat <- function(A.array){
   # Makes an n x pq matrix from a n x p x q array
   q = dim(A.array)[3]
   if (q==1){
      A.mat = F_Sym2Vec(A.array)
   }else{
      A.mat = F_Sym2Vec(A.array[, , 1])
      invisible(sapply(2:q, function(k){A.mat <<- cbind(A.mat, F_Sym2Vec(A.array[, , k]))}))
   }
   return(A.mat)
}
