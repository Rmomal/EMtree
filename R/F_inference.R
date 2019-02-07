
#########################################################################
##############
# test fonction perso
optimRML<-function(init,log.psi, P,eps=1e-6){

  gamma.old=init
  P.vec=F_Sym2Vec(P)
  # while(diff>eps){
  #browser()
  # betaPsi<-F_Vec2Sym(exp(gamma.old))*exp(log.psi)
  # betaPsi[which(betaPsi<1e-16)]=1e-16
  M.vec = F_Sym2Vec(Kirshner(F_Vec2Sym(exp(gamma.old)))$Q)
  lambda = SetLambda(P, F_Vec2Sym(M.vec))
  # browser()
  gamma.vec =log(P.vec/(M.vec+lambda))
  # browser()
  # gamma.vec=gamma.vec-mean(gamma.vec)
  # gamma.vec[which(gamma.vec<(-30))]=-30
  # gamma.vec[which(gamma.vec>30)]=30

  if(sum(is.nan(gamma.vec))!=0) browser()
  gamma.old=gamma.vec

  if(sum(length(which(eigen(Laplacian(F_Vec2Sym(exp(gamma.vec)))[-1,-1])$values<0)))!=0) browser() # le résultat est-il défini positif ?
  # negLik<-(-sum(F_Sym2Vec(P)*(gamma.old+F_Sym2Vec(log.psi)))) + log(SumTree(F_Vec2Sym(exp(gamma.old))))
  #cat("\nneglik= ", negLik," // sum(exp(gamma))",sum(exp(gamma.vec)),"\n")
  #}
  return((gamma.vec))
}

##############
# programme original
F_NegLikelihood <- function(beta.vec, log.psi, P){
  res<-(-sum(F_Sym2Vec(P) * (log(beta.vec) + F_Sym2Vec(log.psi))) + log(SumTree(F_Vec2Sym(beta.vec))))
  return( res)
}
# the derivative of F_NegLikelihood
F_NegGradient <- function(beta.vec, log.psi, P){
  M = Kirshner(F_Vec2Sym(beta.vec)*exp(log.psi))$Q
  lambda = SetLambda(P, M)
  return(- F_Sym2Vec(P)/beta.vec + F_Sym2Vec(M) + rep(lambda, length(beta.vec)))
}
##############
# changement de variable log
F_NegLikelihood_Trans <- function(gamma, log.psi, P){
  if(sum(exp(gamma))>2){
    lambda=1e4
    saveRDS(F_Vec2Sym(exp(gamma)),"beta_negLik.rds")
    print(paste0("Neglik :",sum(exp(gamma))))
  }else{
    M = Kirshner(F_Vec2Sym(exp(gamma))*exp(log.psi))$Q
    lambda = SetLambda(P, M)
    saveRDS(M,"M_negLik.rds")
    saveRDS(F_Vec2Sym(exp(gamma)),"beta_negLik.rds")
    print(paste0("Neglik :",sum(exp(gamma))," // ",sum(P/(M+lambda))," // ", lambda))

  }
  # if(sum(exp(gamma))>2) gamma = log(P/(M+lambda))
  res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+
    log(SumTree(F_Vec2Sym(exp(gamma))))+lambda*(sum(exp(gamma))-0.5)
  if(is.nan(res)){
    browser()
  }
  return( res)
}
F_NegGradient_Trans <- function(gamma, log.psi, P){
  M = Kirshner(F_Vec2Sym(exp(gamma))*exp(log.psi))$Q
  lambda = SetLambda(P, M)
  saveRDS(P,"P.rds")
  saveRDS(M,"M_gradient.rds")
  saveRDS(F_Vec2Sym(exp(gamma)),"beta_gradient.rds")
  print(paste0("Gradient :",sum(exp(gamma))," // ",sum(P/(M+lambda))," // ", lambda))
  # cat( "\nval lambda: ",lambda," sum gamma:",sum(gamma)," sum expGamma:", sum(exp(gamma)))
  return(- F_Sym2Vec(P)+ exp(gamma)*(F_Sym2Vec(M) + lambda))
}
#########################################################################
SetLambda <- function(P, M, eps = 1e-6){
  # F.x has to be increasing. The target value is 0
  #browser()
  F.x <- function(x){
    if(x!=0){
      1 - sum(P / (x+M))
    }else{
      1 - (2*sum(P[upper.tri(P)] / M[upper.tri(M)]))
    }
  }
  x.min = ifelse(F.x(0) >0,-20,1e-4);
  while(F.x(x.min)>0){x.min = x.min -x.min/2}
  x.max = 10
  while(F.x(x.max)<0){x.max = x.max * 2}
  x = (x.max+x.min)/2
  f.min = F.x(x.min)
  f.max = F.x(x.max)
  f = F.x(x)

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
  }
  # print( 1 - sum(P / (x+M)))
  if(abs( 1 - sum(P / (x+M)))>2) browser()
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
FitBetaStatic <- function(beta.init, psi, maxIter, eps1 = 1e-6,eps2=1e-4, optim_method, verbatim,plot){
  options(nwarnings = 1)
  beta.tol = 1e-4
  beta.min = 1e-16
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi)
  iter = 0
  logpY = rep(0, maxIter)
  beta.diff = diff.loglik = 2 * eps2
  cat("\nLikelihoods: ")
  T1<-Sys.time()
  while (((beta.diff > eps1) || (diff.loglik>eps2) ) && iter < maxIter ){
    iter = iter+1

    #print("CALCUL P")
    P = EdgeProba(beta.old*psi)

    init=F_Sym2Vec(beta.old)
    long=length(F_Sym2Vec(beta.old))

    # gamma = optim(log(init), F_NegLikelihood_Trans, gr=F_NegGradient_Trans,method='BFGS', log.psi, P)$par
    #gamma = cma_es(log(init), F_NegLikelihood_Trans,  log.psi, P, lower=rep(1e-30,long), upper=rep(1,long),
    #             control=list(sigma=1e-4))$par

    gamma=optimRML(log(init), log.psi, P, 1e-6)

    beta=exp(gamma)

    beta[which(beta< beta.min)] = beta.min
    beta=F_Vec2Sym(beta)


    diag(beta) =0

    logpY[iter] = -F_NegLikelihood(F_Sym2Vec(beta),log.psi,P)
    cat(logpY[iter],", ")
    beta.diff = max(abs(beta.old-beta))
    beta.old = beta

    if(iter > 1){diff.loglik =  abs(logpY[iter]-logpY[iter-1])}else{diff.loglik=1}

  }

  time<-difftime(Sys.time(),T1)
  logpY = logpY[1:iter]
  g<-tibble(p=logpY) %>% rowid_to_column() %>%
    ggplot(aes(rowid,p))+geom_point()+geom_line()+theme_minimal()+labs(x="Iter",y="Likelihood")
  P = EdgeProba(beta.old*psi)
  if(verbatim){
    cat("\nConvergence took",round(time,2), attr(time, "units")," and ",
        iter," iterations.\nLikelihood difference =", diff.loglik, "\nBetas difference =",beta.diff)
  }
  if(plot) print(g)
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
EMtree<-function(PLNobject,  maxIter, cond.tol=1e-10, optim_method, verbatim=TRUE, plot=FALSE){
  CorY=cov2cor(PLNobject$model_par$Sigma)
  p = ncol(CorY)
  n=PLNobject$n
  alpha.psi = F_AlphaN(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)

  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, maxIter = maxIter, optim_method=optim_method,
                        verbatim=verbatim, plot=plot)

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
                 verbatim=FALSE,plot=FALSE)[c("ProbaCond","maxIter","times")]

    return(inf)
  }, mc.cores=cores)

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
ComparEMtree <- function(counts, vec_covar, O=NULL, v=0.8, B=1e2, maxIter, cond.tol=1e-14,cores=3,f){
  p=ncol(counts)
  split<-strsplit(vec_covar,"\\$")
  models =c("null",split[[1]][2],split[[2]][2],paste(lapply(split, function(x){x[2]}), collapse=" + "))
  cat("\nmodel ",models[1],": \n")
  p1<-ResampleEMtree(counts, "1", O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  cat("\nmodel ",models[2],": \n")
  p2<-ResampleEMtree(counts, vec_covar[1], O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat

  cat("\nmodel ",models[3],": \n")
  p3<-ResampleEMtree(counts, vec_covar[2], O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  cat("\nmodel ",models[4],": \n")
  p4<-ResampleEMtree(counts, vec_covar, O=O, v=v, B=B, maxIter, cond.tol=cond.tol,cores=cores)$Pmat

  Stab.sel=list(p1,p2,p3,p4)

  mat<-data.frame(freq_selec_list(Stab.sel,1,p,f))
  allNets<-tibble(P = list(mat), models =models )  %>%
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
freq_selec<-function(Pmat,p,f){
  return(F_Vec2Sym(1*colMeans(1*(Pmat>2/p))>f))
}

freq_selec_list<-function(list_Pmat,x,p,f){
  return(F_Vec2Sym( 1*(colMeans( 1*(list_Pmat[[x]]>2/p))>f)))
}

select_edges<-function(data,p=p){
  1*(data$ProbaCond>2/p)
}
