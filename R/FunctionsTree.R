######################################################################################
# Summing over trees (Matrix-tree thm)
Laplacian <- function(W){
  # W = squared weight matrix
  return(diag(rowSums(W)) - W)
}

######################################################################################
# Summing over trees (Matrix-tree thm)
SumTree <- function(W){
  # W = squared weight matrix

  return(det(Laplacian(W)[-1,-1]))
}

######################################################################################
# Summing over trees (Matrix-tree thm)
LogSumTree <- function(W){
  # W = squared weight matrix
  return(sum(log(eigen(Laplacian(W)[-1, -1], symmetric=T)$values)))
}

######################################################################################
# Computing edge probability
Kirshner <- function(W){
  # W = squared weight matrix
  # Kirshner (07) formulas
  # W = beta.unif*phi

  if(!isSymmetric(W)){cat('Pb: W non symmetric!')}
  p = nrow(W)
  L = Laplacian(W)[-1, -1]
  if(!is.finite(sum(L))) browser()
  Leigen = eigen(L); Q = (Leigen$vectors) %*% diag(1/Leigen$values) %*% t(Leigen$vectors)
  # Q = chol2inv(chol(L));
  # Q = solve(L[-1, -1]);
  # P = W[-1, -1] * (diag(Q)%o%rep(1, p-1) + rep(1, p-1)%o%diag(Q) - 2*Q)
  # P = rbind(c(0, W[1, -1]*diag(Q)), cbind(W[-1, 1]*diag(Q), P))
  Q = rbind(c(0, diag(Q)),
            cbind(diag(Q), (diag(Q)%o%rep(1, p-1) + rep(1, p-1)%o%diag(Q) - 2*Q)))
  Q = .5*(Q + t(Q))
  P = W * Q
  P = .5*(P + t(P))
  return(list(P=P, Q=Q))
}

KirshnerRM <- function(beta){
  L = Laplacian(beta)[-1, -1]
  # Q<-solve(L)
  Leigen = eigen(L); Q = (Leigen$vectors) %*% diag(1/Leigen$values) %*% t(Leigen$vectors)
  colQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=TRUE)
  rowQ<-matrix(diag(Q),nrow=length(diag(Q)),ncol=length(diag(Q)),byrow=FALSE)
  matrixM<-(colQ+rowQ-2*Q)
  matrixM<-rbind(diag(Q),matrixM)
  matrixM<-cbind(c(0,diag(Q)),matrixM)
  matrixK<-beta*matrixM
  # matrix<-beta*Meila(beta)
  return(list(P=matrixK,Q=matrixM))
}

######################################################################################
# Computing edge probability
EdgeProba <- function(W){
  # W = squared weight matrix
  # Direct calculation
 it=-1
 Wcum = SumTree(W)
  if(!isSymmetric(W)){cat('Pb: W non symmpetric!')}
  while(!is.finite(Wcum)){
    it=it+1
    borne=30-it
    message(cat("W corrected, bound=",borne))

    W.log=log(F_Sym2Vec(W))
    W.center=W.log-mean(W.log)
    W.center[which(W.center<(-borne))]=-borne
    W.center[which(W.center>borne)]=borne
    W=F_Vec2Sym(exp(W.center))
    Wcum = SumTree(W)
  }
  #browser()

  p = nrow(W); P = matrix(0, p, p)

  # cst=1 ;
  # while(min(eigen(W)$values)<0){
  #   W=W+diag(rep(cst,ncol(W)));
  #   cst=1.1*cst
  #  }
  # browser()
  sapply(1:(p-1),
         function(j){
           sapply((j+1):p,
                  function(k){
                    W_jk = W; W_jk[j, k] = W_jk[k, j] = 0 #tuer l'arête kj dans W_kj
                    P[k, j] <<- 1 - SumTree(W_jk) / Wcum
                    P[j, k] <<- P[k, j]
                  }
           )
         }
  )
  if(sum(is.nan(P))!=0) browser()
  P[which(P<1e-10)]=1e-10
  diag(P)=0
  return(P)
}


######################################################################################
# Corrected edge probability, depending on the prior
CorrectedEdgeProba <- function(Pedge, p0=matrix(1/ncol(Pedge), nrow(Pedge), ncol(Pedge)), p=matrix(1/2, nrow(Pedge), ncol(Pedge))){
  return(Pedge*p/p0 / (Pedge*p/p0 + (1-Pedge)*(1-p)/(1-p0)))
}

######################################################################################
# Uniform sampling of a spanning tree
rSpanTree <- function(p){
  A = get.adjacency(mst(graph.adjacency(matrix(runif(p^2), p, p), weighted=TRUE, mode='upper')))
  # A = A + diag(rep(1, p))
  return(as.matrix(A))
}

