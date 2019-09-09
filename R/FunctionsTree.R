
#' Compute a matrix Laplacian
#'
#' @param W Squared weight matrix
#'
#' @return The matrix Laplacian of W
#' @export
#'
#' @examples uniformeW=matrix(1,3,3)
#' Laplacian(uniformeW)
Laplacian <- function(W){
  return(diag(rowSums(W)) - W)
}


#' Summing over trees (Matrix-tree thm)
#'
#' @param W Squared weight matrix
#'
#' @return The number of total possible spanning trees
#' @export
#'
#' @examples uniformeW=matrix(1,3,3)
#' SumTree(uniformeW)
SumTree <- function(W){
  return(det(Laplacian(W)[-1,-1]))
}


#' Computing the M matrix
#'
#' @param W Squared weight matrix
#'
#' @return M matrix as defined in Meila et al. 2006
#' @export
#'
#' @examples W = matrix(c(1,1,3,1,1,1,3,1,1),3,3,byrow=TRUE)
#' Meila(W)
Meila <- function(W){
  if(!isSymmetric(W)){cat('Pb: W non symmetric!')}
  p = nrow(W)
  L = Laplacian(W)[-1, -1]
  Leigen = eigen(L); M = (Leigen$vectors) %*% diag(1/Leigen$values) %*% t(Leigen$vectors)
  M = rbind(c(0, diag(M)),
            cbind(diag(M), (diag(M)%o%rep(1, p-1) + rep(1, p-1)%o%diag(M) - 2*M)))
  M = .5*(M + t(M))
  return(M)
}

#' Computing edges marginal probabilities
#'
#' @param W Squared weight matrix
#'
#' @return Edges probabilities as defined in Kirshner 2007
#' @export
#'
#' @examples W = matrix(c(1,1,3,1,1,1,3,1,1),3,3,byrow=TRUE)
#' Kirshner(W)
Kirshner<-function(W){
  M=Meila(W)
  Pk = W * M
  Pk = .5*(Pk + t(Pk))
  return(Pk)
}


#' Computing edges conditional probabilities
#'
#' @param W squared weight matrix
#'
#' @return Edges conditional probabilities computed directly, without using the Kirshner function.
#' @export
#'
#' @examples W = matrix(c(1,1,3,1,1,1,3,1,1),3,3,byrow=TRUE)
#' EdgeProba(W)
EdgeProba <- function(W){
  it=-1
  Wcum = SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  while(!is.finite(Wcum)){
    #handles numerical issues with matrix tree theorem
    it=it+1
    borne=30-it
    message(cat("W corrected, bound=",borne))

    W.log=log(F_Sym2Vec(W))
    W.center=W.log-mean(W.log)
    W.center[which(W.center<(-borne))]=-borne
    W=F_Vec2Sym(exp(W.center))
    Wcum = SumTree(W)
  }

  p = nrow(W); P = matrix(0, p, p)
  #core of computation
  sapply(1:(p-1),
         function(j){
           sapply((j+1):p,
                  function(k){
                    W_jk = W; W_jk[j, k] = W_jk[k, j] = 0 #kills kj edge in W_kj
                    P[k, j] <<- 1 - SumTree(W_jk) / Wcum
                    P[j, k] <<- P[k, j]
                  }
           )
         }
  )
  P[which(P<1e-10)]=1e-10
  diag(P)=0
  return(P)
}
