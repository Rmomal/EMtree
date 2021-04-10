
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


#' Summing over trees (Matrix Tree Theorem)
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
  p <- nrow(W) ; index<-1
  L <- Laplacian(W)[-index,-index]
  Mei <-solve(L)
  Mei <- rbind(c(0, diag(Mei)),
              cbind(diag(Mei),
                    (diag(Mei) %o% rep(1, p - 1) +
                       rep(1, p - 1) %o% diag(Mei) - 2 * Mei)
              )
  )
  Mei<- 0.5 * (Mei + t(Mei))
  return(Mei)
}

#' Computing edge probabilities using Kirshner (07) formula
#'
#' @param W Squared weight matrix
#'
#' @return Edges probabilities as defined in Kirshner 2007
#' @export
#' @examples W = matrix(c(1,1,3,1,1,1,3,1,1),3,3,byrow=TRUE)
#' Kirshner(W)
Kirshner<-function(W){
  # W = squared weight matrix
  # Kirshner (07) formulas
  p <- nrow(W)
  L <- Laplacian(W)[-1,-1]
  #no need for gmp improved capacity thanks to the adaptive
  # conditioning of psi and beta
  Q<-solve(L)
  Q <- rbind(c(0, diag(Q)),
            cbind(diag(Q), (diag(Q)%o%rep(1, p-1) +
                              rep(1, p-1)%o%diag(Q) - 2*Q)))
  Q <- .5*(Q + t(Q))
  P <- W * Q
  P <- .5*(P + t(P))
  return(P)
}


#' Computing edge probabilities
#'
#' @param W squared weight matrix
#' @param verbatim controls verbosity
#'
#' @return Edges conditional probabilities computed directly,
#' without using the Kirshner function.
#' @export
#'
#' @examples W = matrix(c(1,1,3,1,1,1,3,1,1),3,3,byrow=TRUE)
#' EdgeProba(W)
EdgeProba <- function(W, verbatim=FALSE){
  it<--1
  Wcum <- SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  while(!is.finite(Wcum)){
    #handles numerical issues with matrix tree theorem
    it<-it+1
    borne<-30-it
   if(verbatim) message(cat("W corrected, bound=",borne))
    W.log<-log(ToVec(W))
    W.center<-W.log-mean(W.log)
    W.center[which(W.center<(-borne))]<--borne
    W<-ToSym(exp(W.center))
    Wcum <- SumTree(W)
  }

  p <- nrow(W); P <- matrix(0, p, p)
  #core of computation
  invisible(lapply(1:(p-1),
         function(j){
          invisible(lapply((j+1):p,
                  function(k){#kills kj edge in W_kj
                    W_jk <- W; W_jk[j, k] <- W_jk[k, j] <- 0
                    P[k, j] <<- 1 - SumTree(W_jk) / Wcum
                    P[j, k] <<- P[k, j]
                  }
           ))
         }
  ))
  P[which(P<1e-10)]<-1e-10
  diag(P)<-0
  return(P)
}
