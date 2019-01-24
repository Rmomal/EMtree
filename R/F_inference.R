
#########################################################################

F_NegLikelihood <- function(beta.vec, log.psi, P){
  options(warn=-1)
  return( suppressWarnings(-sum(F_Sym2Vec(P) * (log(beta.vec) + F_Sym2Vec(log.psi))) + log(SumTree(F_Vec2Sym(beta.vec)))))
}

# the derivative of F_NegLikelihood
F_NegGradient <- function(beta.vec, log.psi, P){
  M = Kirshner(F_Vec2Sym(beta.vec)*exp(log.psi))$Q
  lambda = SetLambda(P, M)
  return(- F_Sym2Vec(P)/beta.vec + F_Sym2Vec(M) + rep(lambda, length(beta.vec)))
}
#########################################################################
SetLambda <- function(P, M, eps = 1e-6){
  # F.x has to be increasing. The target value is 0
  F.x <- function(x){
    if(x!=0){
      1 - sum(P / (x+M))
    }else{
      1 - (2*sum(P[upper.tri(P)] / M[upper.tri(M)]))
    }
  }
  # suite=TRUE
  # if(F.x(1e-16) >0){
  #   F.x <- function(x){0.99 - sum(P / (x+M))}
  #   if(F.x(1e-16) >0){
  #     suite=FALSE
  #     x=NA
  #     browser()
  #     cat("ECHEC GRADIENT DESCENT : F.X(1e-16) positive at value",F.x(1e-16)," \n")
  #   }
  # }
  # if(suite){
  x.min = ifelse(F.x(0) >0,-20,1e-4);
  while(F.x(x.min)>0){x.min = x.min -x.min/2}
  x.max = 10
  while(F.x(x.max)<0){x.max = x.max * 2}
  x = (x.max+x.min)/2


  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)
  # x.list = exp(seq(log(x.min), log(x.max), length.out=50))
  # plot(x.list, sapply(x.list, function(l){F.x(l)}), type='l', log='xy'); abline(h=0)
  # points(c(x.min, x, x.max), c(f.min, f, f.max), col=c(1, 2, 1))
  # browser()
  while(abs(x.max-x.min) > eps){
    if(f > 0) {
      x.max = x
      f.max = f
    } else{
      x.min = x
      f.min = f
    }
    x = (x.max+x.min)/2;
    f = F.x(x)

    # points(c(x.min, x, x.max), c(f.min, f, f.max), col=c(1, 2, 1))
    # cat(x.min, x, x.max, '/', f.min, f, f.max, '\n')
  }
  # cat("lambda : ", x,"F.x(lambda) : ",f,"\n")
  #x }
  return(x)
}


#########################################################################
# Choice of alpha for alpha * n
F_AlphaN <- function(CorY, n, cond.tol=1e-10){
  # Grid on alpha
  alpha.grid = (1:n)/n; alpha.nb = length(alpha.grid);
  cond = Inf; a = 0
  while(cond > cond.tol && a<length(alpha.grid)){
    a = a+1
    psi.vec = F_Sym2Vec(-alpha.grid[a]*n*log(1 - CorY^2)/2);
    psi.vec = psi.vec - mean(psi.vec)
    psi = F_Vec2Sym(exp(psi.vec))
    lambda = svd(psi)$d
    cond = min(abs(lambda))/max(abs(lambda))
  }
  alpha = alpha.grid[a-1]
  psi.vec = F_Sym2Vec(-alpha*n*log(1 - CorY^2)/2);
  psi.vec = psi.vec - mean(psi.vec)
  psi = F_Vec2Sym(exp(psi.vec))
  return(list(psi=psi, alpha=alpha))
}

#########################################################################
FitBetaStatic <- function(beta.init, psi, maxIter, eps1 = 1e-6,eps2=1e-4, optim_method, verbatim){
  options(nwarnings = 1)
  beta.tol = 1e-4
  beta.min = 1e-30
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi)
  iter = 0
  logpY = rep(0, maxIter)
  beta.diff = diff.loglik = 2 * eps2

  T1<-Sys.time()
  while (((beta.diff > eps1) || (diff.loglik>eps2) ) && iter < maxIter ){
    iter = iter+1

    P = EdgeProba(beta.old*psi)
    beta =  tryCatch(F_Vec2Sym(
      optim(F_Sym2Vec(beta.old), F_NegLikelihood, gr=F_NegGradient,method="BFGS",
            log.psi, P)$par),
      warning = function(w) {print("Negative weights during optimization process" )}
    )

    #browser()
    beta[which(beta< beta.min)] = beta.min
    diag(beta) = 0
    logpY[iter] = -F_NegLikelihood(F_Sym2Vec(beta),log.psi,P)

    beta.diff = max(abs(beta.old-beta))
    beta.old = beta

    if(iter > 1){diff.loglik =  abs(logpY[iter]-logpY[iter-1])}else{diff.loglik=1}

  }
  time<-difftime(Sys.time(),T1)
  logpY = logpY[1:iter]
  P = EdgeProba(beta.old*psi)
  if(verbatim){
    print(paste0("Convergence took",round(time,2), attr(time, "units")," and ",
                 iter," iterations.\nLikelihood difference =", diff.loglik, "\nBetas difference =",beta.diff))
  }
  return(list(beta=beta, logpY=logpY,ProbaCond=P,maxIter=iter, times=time))
}

#########################################################################
#' Title
#'
#' @param PLNobject
#' @param maxIter
#' @param cond.tol
#' @param optim_method
#' @param verbatim
#'
#' @return
#' @export
#'
#' @examples
EMtree<-function(PLNobject,  maxIter, cond.tol=1e-10, optim_method, verbatim=TRUE){
  CorY=cov2cor(PLNobject$model_par$Sigma)
  p = ncol(CorY)
  n=PLNobject$n
  alpha.psi = F_AlphaN(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)

  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, maxIter = maxIter, optim_method=optim_method,
                        verbatim=verbatim)

  return(FitEM)
}

#################################################################################
# Resampling edge probability
# Y, X, O: same as for PLN
# v = (subsample size) / (total sample size)
# B = nb resamples
# Out = Pmat = B x p(p-1)/2 matrix with edge probability for each resample
# Y = Data$count; X = matrix(1, n, 1); O = matrix(0, n, p)
#' Title
#'
#' @param counts
#' @param vec_covar
#' @param O
#' @param v
#' @param B
#' @param maxIter
#' @param cond.tol
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
ResampleEMtree <- function(counts, vec_covar, O=NULL, v=0.8, B=1e2, maxIter, cond.tol=1e-14,cores=3){
  n = nrow(counts)
  p = ncol(counts)
  P = p * (p - 1) / 2
  V = round(v * n)
  Pmat = matrix(0, B, P)
  if(is.null(O)){ O=matrix(0, n, p)}

  string<-paste("counts", paste(vec_covar, collapse=" + "), sep=" ~ ")
  formula<-as.formula(string)
  X = as.matrix(lm(formula, x=T)$x)

  obj<-mclapply(1:B,function(b){
    set.seed(b)
    sample = sample(1:n, V, replace = F)
    counts.sample = counts[sample,]
    X.sample = X[sample,]
    O.sample = O[sample,]

    suppressWarnings(
      PLN.sample <- PLN(counts.sample ~ -1 + X.sample + offset(O.sample),control = list("trace"=0))
    )

    inf<-EMtree( PLN.sample, maxIter=maxIter, cond.tol=cond.tol,
                 verbatim=FALSE)[c("ProbaCond","maxIter","times")]

    return(inf)
  },mc.cores=cores)

  Pmat<-do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$ProbaCond)}))
  summaryiter = do.call(c,lapply(obj,function(x){x$maxIter}))
  times<-do.call(c,lapply(obj,function(x){x$time}))

  return(list(Pmat=Pmat,maxIter=summaryiter,times=times))
}

#' Title
#'
#' @param counts
#' @param vec_covar
#' @param O
#' @param v
#' @param B
#' @param maxIter
#' @param cond.tol
#' @param cores
#' @param f
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
ComparEMtree <- function(counts, vec_covar, O=NULL, v=0.8, B=1e2, maxIter, cond.tol=1e-14,cores=3,f,seed){
  p=ncol(counts)
  p1<-ResampleEMtree(counts, "1", O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  p2<-ResampleEMtree(counts, vec_covar[1], O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  p3<-ResampleEMtree(counts, vec_covar[2], O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  p4<-ResampleEMtree(counts, vec_covar, O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  Stab.sel=list(p1,p2,p3,p4)

  mat<-data.frame(freq_selec_list(Stab.sel,1,p,f))
  set.seed(seed)
  allNets<-tibble(P = list(mat), models =c("null",vec_covar[1],vec_covar[2],paste(vec_covar, collapse=" + ")) )  %>%
    mutate(P=map( seq_along(P), function(x) {
      df<-freq_selec_list(Stab.sel,x,p,f)
      df[lower.tri(df, diag = TRUE)]<-0
      df<-data.frame(df)
      colnames(df)<-1:ncol(df)
      df
    })) %>%
    mutate(P = map(P,~rownames_to_column(.) %>%
                     gather(key, value , -rowname))) %>%
    unnest()
  allNets<-allNets[,c(3,2,1,4)]
  return(allNets)
}


#' Title
#'
#' @param list
#' @param p
#' @param f
#'
#' @return
#' @export
#'
#' @examples
freq_selec<-function(list,p,f){
 return(F_Vec2Sym(1*colMeans(1*(list>2/p))>f))
}




freq_selec_list<-function(list,x,p,f){
  return(F_Vec2Sym( 1*(colMeans( 1*(list[[x]]>2/p))>f)))
}
select_edges<-function(data,p=p){
  1*(data$ProbaCond>2/p)
}
