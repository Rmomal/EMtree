
#########################################################################
##############
source("R/FunctionsMatVec.R")
source("R/FunctionsTree.R")

##############
# gradients with log change of variable
F_NegLikelihood_Trans <- function(gamma, log.psi, P){
  gamma=gamma-mean(gamma)
  gamma[which(gamma<(-30))]=-30
  M = Kirshner(F_Vec2Sym(exp(gamma)))$Q
  lambda = SetLambda(P, M)

  res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+log(SumTree(F_Vec2Sym(exp(gamma))))+
    lambda*(sum(exp(gamma))-0.5)
  if(is.nan(res)){
    cat(max(gamma),": higher bound ")
    gamma[which(gamma>(30))]=30
    M = Kirshner(F_Vec2Sym(exp(gamma)))$Q
    lambda = SetLambda(P, M)

    res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+log(SumTree(F_Vec2Sym(exp(gamma))))+
      lambda*(sum(exp(gamma))-0.5)

    if(is.nan(res)) browser()
  }

  return( res)
}
F_NegGradient_Trans <- function(gamma, log.psi, P){
  M = Kirshner(F_Vec2Sym(exp(gamma)))$Q
  lambda = SetLambda(P, M)
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
FitBetaStatic <- function(beta.init, psi, maxIter, eps1 = 1e-6,eps2=1e-4, verbatim,plot){
  options(nwarnings = 1)
  beta.tol = 1e-4
  beta.min = 1e-16
  beta.old = beta.init / sum(beta.init)
  log.psi = log(psi)
  iter = 0
  logpY = rep(0, maxIter)
  beta.diff = diff.loglik = 2 * eps2
  if(verbatim) cat("\nLikelihoods: ")
  T1<-Sys.time()
  stop=FALSE
  while (((beta.diff > eps1) || (diff.loglik>eps2) ) && iter < maxIter && !stop){
    iter = iter+1
    P = EdgeProba(beta.old*psi)
    init=F_Sym2Vec(beta.old)
    long=length(F_Sym2Vec(beta.old))
    browser()
    gamma = optim(log(init), F_NegLikelihood_Trans, gr=F_NegGradient_Trans,method='BFGS', log.psi, P)$par
    beta=exp(gamma)
    beta[which(beta< beta.min)] = beta.min
    beta=F_Vec2Sym(beta)

    logpY[iter] = -F_NegLikelihood(F_Sym2Vec(beta),log.psi,P)
    if(verbatim) cat(logpY[iter],", ")
    beta.diff = max(abs(beta.old-beta))
    beta.old = beta
    diffPres=logpY[iter]-logpY[iter-1]
    if(iter > 1){diff.loglik =  abs(diffPres)}else{diff.loglik=1}
  }

  time<-difftime(Sys.time(),T1)
  logpY = logpY[1:iter]
  g<-tibble(p=logpY) %>% rowid_to_column() %>%
    ggplot(aes(rowid,p))+geom_point()+geom_line()+theme_minimal()+labs(x="Iter",y="Likelihood")
  P = EdgeProba(beta.old*psi)
  if(verbatim){
    cat("\nConvergence took",round(time,2), attr(time, "units")," and ",
        iter," iterations.\nLikelihood difference =", diff.loglik, "\nBetas difference =",beta.diff)
  }else{
    cat("\nConvergence took",round(time,2), attr(time, "units")," and ",
        iter," iterations.")
  }
  if(plot) print(g)
  return(list(beta=beta, logpY=logpY,ProbaCond=P,maxIter=iter, timeEM=time))
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
EMtree<-function(PLNobject,  maxIter=20, cond.tol=1e-10, verbatim=TRUE, plot=FALSE){
  CorY=cov2cor(PLNobject$model_par$Sigma)
  p = ncol(CorY)
  n=PLNobject$n
  alpha.psi = F_AlphaN(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  print(alpha.psi$alpha)
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)

  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, maxIter = maxIter,
                        verbatim=verbatim, plot=plot)
  FitEM$alpha=alpha.psi$alpha
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
#' @param counts Data of observed counts with dimensions n x p, either a matrix, data.frame or tibble.
#' @param vec_covar Vecteur of quoted covariates to be accounted for with dimension n x d, with the data.frame of covariates. For example with
#' X matrix of covariates : c("X$first", "X$second").
#' @param O The n x p matrix of offsets, filled with zeros by default.
#' @param v The proportion of observed data to be taken in each sub-sample. It is the ratio (sub-sample size)/n
#' @param S Total number of sub-samples.
#' @param maxIter Maximum number of EMtree iterations at each sub-sampling.
#' @param cond.tol Tolerance for the psi matrix.
#' @param cores Number of cores, can be greater than 1 if data involves less than about 32 species.
#'
#' @return Returns a list which contains the Pmat data.frame, and vectors of EMtree maximum iterations and running times in each
#' resampling.
#'\itemize{
#'  \item{Pmat: }{S x p(p-1)/2 matrix with edge probabilities for each resample}
#'  \item{maxIter: }{EMtree maximum iterations in each resampling.}
#'  \item{times: }{EMtree running times in each resampling.}
#' }
#'
#' @export
#'
#' @examples
ResampleEMtree <- function(counts, vec_covar=NULL,data_covar=NULL, covariate=NULL  , O1=NULL, O2=NULL, v=0.8, S=1e2, maxIter, cond.tol=1e-14,cores=3){
  #browser()
  counts=as.matrix(counts)
  n = nrow(counts)
  p = ncol(counts)
  P = p * (p - 1) / 2
  V = round(v * n)
  Pmat = matrix(0, S, P)
  if(is.null(O1)){ O1=matrix(1, n, p)}
  if(is.null(O2)){ O2=matrix(1, n, p)}
  #browser()
  if(is.null(covariate)){
    attach(data_covar)
    string<-paste("counts", paste(vec_covar, collapse=" + "), sep=" ~ ")
    formula<-as.formula(string)

    X = as.matrix(lm(formula, x=T)$x)
  }else{
    X=covariate
  }

  obj<-mclapply(1:S,function(b){
    cat("\nS=",b," ")
    set.seed(b)

    sample = sample(1:n, V, replace = F)
    counts.sample = counts[sample,]
    X.sample = X[sample,]
    O1.sample = O1[sample,]
    O2.sample = O2[sample,]

    suppressWarnings(
      PLN.sample <- PLN(counts.sample ~ -1 + X.sample + offset(log(O1.sample))+ offset(log(O2.sample)),control = list("trace"=0))
    )

    inf1<-EMtree( PLN.sample, maxIter=maxIter, cond.tol=cond.tol,
                  verbatim=FALSE,plot=FALSE)[c("ProbaCond","maxIter","timeEM","alpha")]
    cat(" ",inf1$alpha)

    inf<-inf1[c("ProbaCond","maxIter","timeEM")]
    return(inf)
  }, mc.cores=cores)

  lengths<- sapply(obj,function(x){length(x)})
  if(mean(lengths)!=3){
    indices<-which(lengths!=3)
    lapply(indices, function(x){
      set.seed(x)
      sample = sample(1:n, V, replace = F)
      counts.sample = counts[sample,]
      X.sample = X[sample,]
      O.sample = O[sample,]

      suppressWarnings(
        PLN.sample <- PLN(counts.sample ~ -1 + X.sample + offset(O.sample),control = list("trace"=0))
      )
      obj[[x]]<<-EMtree( PLN.sample, maxIter=maxIter, cond.tol=cond.tol,
                         verbatim=FALSE,plot=FALSE)[c("ProbaCond","maxIter","timeEM")]
    })

  }
  Pmat<-do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$ProbaCond)}))

  summaryiter = do.call(c,lapply(obj,function(x){x$maxIter}))
  times<-do.call(c,lapply(obj,function(x){x$timeEM}))

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
ComparEMtree <- function(counts, vec_covar, O=NULL, v=0.8, S=1e2, maxIter, cond.tol=1e-14,cores=3,f){
  p=ncol(counts)
  # browser()
  split<-strsplit(vec_covar,"\\$")
  models =c("null",split[[1]][2],split[[2]][2],paste(lapply(split, function(x){x[2]}), collapse=" + "))
  cat("\nmodel ",models[1],": \n")
  p1<-ResampleEMtree(counts, vec_covar="1",covariate=NULL, O=O, v=v, S=S, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  cat("\nmodel ",models[2],": \n")
  p2<-ResampleEMtree(counts, vec_covar=vec_covar[1],covariate=NULL, O=O, v=v, S=S, maxIter, cond.tol=cond.tol,cores=cores)$Pmat

  cat("\nmodel ",models[3],": \n")
  p3<-ResampleEMtree(counts, vec_covar=vec_covar[2], covariate=NULL,O=O, v=v, S=S, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  cat("\nmodel ",models[4],": \n")
  p4<-ResampleEMtree(counts, vec_covar=vec_covar,covariate=NULL, O=O, v=v, S=S, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  #browser()
  Stab.sel=list(p1,p2,p3,p4)

  mat<-data.frame(freq_selec_list(Stab.sel,1,p,f)) # the first element is initialized
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
  #  browser()
  return(F_Vec2Sym( 1*(colMeans( 1*(list_Pmat[[x]]>2/p))>f)))
}

select_edges<-function(data,p=p){
  1*(data$ProbaCond>2/p)
}

freq_selec_pmat<-function(list,p,f){ # rend une liste de sélection pour chaque elt (modeles) de la liste de départ
  lapply(list, function(x){
    1*(colMeans(1*(x$Pmat>2/p))>f)
  })
}
