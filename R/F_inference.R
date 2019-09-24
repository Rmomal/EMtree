
#########################################################################
##############
source("R/FunctionsMatVec.R")
source("R/FunctionsTree.R")

##############
F_NegLikelihood <- function(beta.vec, log.psi, P){
  M = Meila(F_Vec2Sym(beta.vec))
  lambda = SetLambda(P, M)
  return(- sum(F_Sym2Vec(P)*(log(beta.vec)+F_Sym2Vec(log.psi))) +
           log(SumTree(F_Vec2Sym(beta.vec)))+
           lambda*(sum(beta.vec)-0.5))
}
# gradients with log change of variable
F_NegLikelihood_Trans <- function(gamma, log.psi, P, trim=TRUE){
  if(trim){
    gamma=gamma-mean(gamma)
    gamma[which(gamma<(-30))]=-30
  }
  M = Meila(F_Vec2Sym(exp(gamma)))
  lambda = SetLambda(P, M)

  res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+
    log(SumTree(F_Vec2Sym(exp(gamma))))+
    lambda*(sum(exp(gamma))-0.5)

  if(is.nan(res)){
    cat(max(gamma),": higher bound ")
    gamma[which(gamma>(30))]=30
    M = Meila(F_Vec2Sym(exp(gamma)))
    lambda = SetLambda(P, M)
    res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+
      log(SumTree(F_Vec2Sym(exp(gamma))))+
      lambda*(sum(exp(gamma))-0.5)
    if(is.nan(res)) browser()
  }
  return( res)
}

F_NegGradient_Trans <- function(gamma, log.psi, P){
  M = Meila(F_Vec2Sym(exp(gamma)))
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
#
#' Choice of alpha for a numerically stable psi matrix
#'
#' @param CorY Correlation matrix of dimensoins p x p with p the number of species for example.
#' @param n number of samples
#' @param cond.tol tolerance for the matrix conditionment, which is here the ratio of the absolute values of minimal eigen value on maximal eigen value
#'
#' @return \itemize{
#' \item{psi: }{p x p corrected psi matrix}
#' \item{alpha: }{final alpha value needed to correct the psi matrix}}
#' @export
#'
#' @examples n=30
#' p=10
#' S=5
#' Y=data_from_scratch("tree",p=p,n=n)$data
#' beta = matrix(1/10,10,10)
#' psi=Psi_alpha(cor(Y), n)$psi
Psi_alpha <- function(CorY, n, cond.tol=1e-10){
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
#' internal function
#'
#' @param beta.init initialization of beta weights
#' @param psi psi matrix, filled with ratios of bivariate probabilities over marginals, which can in the Gaussian
#' wase be deduced from the correlation matrix.
#' @param maxIter maximum number of iterations
#' @param eps1 variation higher bound of beta weights
#' @param eps2 variation higher bound of log likelihood
#' @param verbatim boolean
#' @param plot boolean
#'
#' @return
#' \itemize{
#'  \item{edges_prob: }{p x p matrix of edges probabilities}
#'  \item{edges_weight: }{p x p matrix of edges weights for any spanning tree}
#'  \item{logpY: }{vector of log-likelihoods}
#'  \item{maxIter: }{final number of iterations EMtree has ran}
#'  \item{timeEM: }{EMtree computation time}
#' }
#' @export
#' @importFrom tibble tibble rowid_to_column
#' @importFrom ggplot2 ggplot geom_point geom_line theme_minimal labs aes
#' @examples n=30
#' p=10
#' S=5
#' Y=data_from_scratch("tree",p=p,n=n)$data
#' beta = matrix(1/10,10,10)
#' psi=Psi_alpha(cor(Y), n)$psi
#' FitEM = FitBetaStatic(beta.init=beta, psi=psi, maxIter = 6, verbatim=TRUE, plot=TRUE)
FitBetaStatic <- function(beta.init, psi, maxIter=50, eps1 = 1e-6,eps2=1e-4, verbatim=TRUE,plot=FALSE){
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
  if(plot){

    g<-tibble(p=logpY) %>% rowid_to_column() %>%
      ggplot(aes(rowid,p))+geom_point()+geom_line()+theme_minimal()+labs(x="Iter",y="Likelihood")
    print(g)
  }
  P = EdgeProba(beta.old*psi)
  if(verbatim){
    cat("\nConvergence took",round(time,2), attr(time, "units")," and ",
        iter," iterations.\nLikelihood difference =", diff.loglik, "\nBetas difference =",beta.diff)
  }else{
    cat("\nConvergence took",round(time,2), attr(time, "units")," and ",
        iter," iterations.")
  }

  return(list(edges_prob=P, edges_weight=beta, logpY=logpY,maxIter=iter, timeEM=time))
}

#########################################################################
#' Computes edges probability
#'
#' @param PLNobject  an object resulting from the use of the `PLN` function from package `PLNmodels`
#' @param maxIter Maximum number of iterations for EMtree
#' @param cond.tol Tolerence parameter for the conditionement of psi matrix
#' @param verbatim talks if set to TRUE
#' @param plot plots likelihood if set to TRUE
#'
#' @return
#' \itemize{
#'  \item{edges_prob: }{p x p matrix of edges probabilities}
#'  \item{edges_weight: }{p x p matrix of edges weights for any spanning tree}
#'  \item{logpY: }{vector of log-likelihoods}
#'  \item{maxIter: }{final number of iterations EMtree has ran}
#'  \item{timeEM: }{EMtree computation time}
#'  \item{alpha: }{data signal/noise ratio}
#' }
#' @export
#' @examples
#' n=30
#' p=10
#' Y=data_from_scratch("tree",p=p)$data
#' PLN_Y = PLNmodels::PLN(Y~1)
#' EMtree(PLN_Y,verbatim=TRUE, plot=TRUE)
EMtree<-function(PLNobject,  maxIter=30, cond.tol=1e-10, verbatim=TRUE, plot=FALSE){
  CorY=cov2cor(PLNobject$model_par$Sigma)
  p = ncol(CorY)
  n=PLNobject$n
  alpha.psi = Psi_alpha(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi
  print(alpha.psi$alpha)
  beta.unif = matrix(1, p, p); diag(beta.unif) = 0; beta.unif = beta.unif / sum(beta.unif)

  FitEM = FitBetaStatic(beta.init=beta.unif, psi=psi, maxIter = maxIter,
                        verbatim=verbatim, plot=plot)
  FitEM$alpha=alpha.psi$alpha
  return(FitEM)
}

#################################################################################
#' Resampling procedure for  edges probability
#'
#' @param counts Data of observed counts with dimensions n x p, either a matrix, data.frame or tibble.
#' @param covar_matrix matrix of covariates, should have the same number of rows as the count matrix.
#' @param O Matrix of offsets, with dimension n x p
#' @param v The proportion of observed data to be taken in each sub-sample. It is the ratio (sub-sample size)/n
#' @param S Total number of wanted sub-samples.
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
#' @importFrom PLNmodels PLN
#' @importFrom parallel mclapply
#' @examples
#'n=30
#'p=10
#'S=5
#'Y=data_from_scratch("tree",p=p,n=n)$data
#'X = data.frame(rnorm(n),rbinom(n,1,0.7))
#'resample=ResampleEMtree(Y,covar_matrix=X, S=S,cores = 1)
#'str(resample)
ResampleEMtree <- function(counts, covar_matrix=NULL  , O=NULL, v=0.8, S=1e2, maxIter=30, cond.tol=1e-14,cores=3){

  counts=as.matrix(counts)
  n = nrow(counts)
  p = ncol(counts)
  P = p * (p - 1) / 2
  V = round(v * n)
  Pmat = matrix(0, S, P)
  if(is.null(O)){ O=matrix(1, n, p)}

  if(is.null(covar_matrix)){
    #default
    X=matrix(1,nrow=n,ncol=1)
  }else{
    X=as.matrix(covar_matrix)
  }

  obj<-mclapply(1:S,function(b){
    cat("\nS=",b," ")
    set.seed(b)

    sample = sample(1:n, V, replace = F)
    counts.sample = counts[sample,]
    X.sample = data.frame(X[sample,])
    O.sample = O[sample,]

    suppressWarnings(
      PLN.sample <- PLN(counts.sample ~ -1  + offset(log(O.sample)) + ., data=X.sample, control = list("trace"=0))
    )

    inf1<-EMtree( PLN.sample, maxIter=maxIter, cond.tol=cond.tol,
                  verbatim=FALSE,plot=FALSE)[c("edges_prob","maxIter","timeEM","alpha")]
    cat(" ",inf1$alpha)

    inf<-inf1[c("edges_prob","maxIter","timeEM")]
    return(inf)
  }, mc.cores=cores)
  #parallelization can malfunction.
  #here chack if all results were correctly computed (correct number of results)
  #if not compute the missing results sequentially
  lengths<- sapply(obj,function(x){length(x)})
  if(mean(lengths)!=3){
    indices<-which(lengths!=3)
    lapply(indices, function(x){
      set.seed(x)
      sample = sample(1:n, V, replace = F)
      counts.sample = counts[sample,]
      X.sample = data.frame(X[sample,])
      O.sample = O[sample,]

      suppressWarnings(
        PLN.sample <- PLN(counts.sample ~ -1  + offset(log(O.sample)) + ., data=X.sample, control = list("trace"=0))
      )
      obj[[x]]<<-EMtree( PLN.sample, maxIter=maxIter, cond.tol=cond.tol,
                         verbatim=FALSE,plot=FALSE)[c("edges_prob","maxIter","timeEM")]
    })

  }

  Pmat<-do.call(rbind,lapply(obj,function(x){F_Sym2Vec(x$edges_prob)}))

  summaryiter = do.call(c,lapply(obj,function(x){x$maxIter}))
  times<-do.call(c,lapply(obj,function(x){x$timeEM}))

  return(list(Pmat=Pmat,maxIter=summaryiter,times=times))
}


#' Runs EMtree for several covariates choices
#'
#' @param counts Data of observed counts with dimensions n x p, either a matrix, data.frame or tibble.
#' @param covar_matrix matrix of covariates, should have the same number of rows as the count matrix.
#' @param models list of covariate combinations to be tested. For example list(1,2) will design two linear models with the first two covariates adjusted separately
#' @param m_names list of names for the models to be compared, for example list("first model","last model")
#' @param O Offset matrix with dimensions n x p
#' @param Pt Probability threshold for every sub-sample
#' @param v The proportion of observed data to be taken in each sub-sample. It is the ratio (sub-sample size)/n
#' @param S Number of desired sub-samples
#' @param maxIter Maximum number of iterations of EMtree
#' @param cond.tol Tolerance for the psi matrix.
#' @param cores Number of cores, can be greater than 1 if data involves less than about 32 species.
#'
#' @return a tibble in a suitable format for graphical purposes. It has four columns describing each edge of each model: the origin node ("node1"),
#' the destination node ("node2"), the model it corresponds to ("model") and its weight ("value").
#' @export
#' @importFrom parallel mclapply
#' @importFrom purrr map
#' @importFrom tidyr gather  unnest
#' @importFrom dplyr mutate filter
#' @importFrom tibble tibble rownames_to_column
#' @examples
#'n=30
#'p=10
#'S=3
#'Y = data_from_scratch(type="tree",p=p,n=n)$data
#'X = data.frame(rnorm(n),rbinom(n,1,0.7))
#'ComparEMtree(Y,X,models=list(1,2,c(1,2)),m_names=list("1","2","both"),Pt=0.3,S=S, cores=1)
ComparEMtree <- function(counts, covar_matrix, models, m_names, O=NULL, Pt=0.1, v=0.8, S=1e2, maxIter=50, cond.tol=1e-14,cores=3){

  Stab.sel<- lapply(seq_along(models),function(x){
    cat("\nmodel ",m_names[[x]])
    ResampleEMtree(counts,covar_matrix[,models[[x]]],O=O, v=v, S=S, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
  })
  p=ncol(counts)

  mat<-data.frame(freq_selec(Stab.sel[[1]],Pt)) # the first element is initialized
  allNets<-tibble(P = list(mat), mods =m_names )  %>%
    mutate(P=map( seq_along(P), function(x) {
      df<-freq_selec(Stab.sel[[x]],Pt=Pt)
      df[lower.tri(df, diag = TRUE)]<-NA
      df<-data.frame(df)
      colnames(df)<-1:ncol(df)
      df
    })) %>%
    mutate(P = map(P,~rownames_to_column(.) %>%
                     gather(key, value , -rowname) %>%filter(!is.na(value))
                   ),
           mods=unlist(mods)
    ) %>%
    unnest(cols = c(P))
  allNets<-allNets[,c(1,2,4,3)]
  colnames(allNets) = c("node1","node2","model","weight")
  return(allNets)
}

####################################################


#' Computes edges selection frequency after resampling procedure
#'
#' @param Pmat matrix gathering edgges probability, with dimensions number of resamples x number of possible edges. Typically the Pmat output of ResampleEMtree()
#' @param Pt edges probabilty threshold
#'
#' @return p x p matrix of edges selection frequency
#' @export
#'
#' @examples
#'n=30
#'p=10
#'S=5
#'Y=data_from_scratch("tree",p=p)$data
#'resample_prob=ResampleEMtree(Y, S=S,cores = 1)$Pmat
#'edges_freq<-freq_selec(resample_prob,Pt=0.3)
#'str(edges_freq)

freq_selec<-function(Pmat,Pt){
  E=ncol(Pmat)
  p=sqrt(2*E+(1/4))-1/2
  if(is.null(Pt)) Pt= 2/p + 0.1
  return(F_Vec2Sym(colMeans(1*(Pmat>Pt))))
}

select_edges<-function(data,p=p){
  1*(data$ProbaCond>2/p)
}

freq_selec_pmat<-function(list,p,f){ # rend une liste de sélection pour chaque elt (modeles) de la liste de départ
  lapply(list, function(x){
    1*(colMeans(1*(x$Pmat>2/p))>f)
  })
}
