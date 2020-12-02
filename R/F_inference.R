
F_NegLikelihood <- function(beta.vec, log.psi, P,sum.constraint){
  M = Meila(F_Vec2Sym(beta.vec))
  lambda = SetLambda(P, M,sum.constraint)
  return(- sum(F_Sym2Vec(P)*(log(beta.vec+(beta.vec==0))+F_Sym2Vec(log.psi))) +
           log(SumTree(F_Vec2Sym(beta.vec)))+
           lambda*(sum(beta.vec)-sum.constraint/2))
}

F_NegGradient_Trans <- function(gamma, log.psi, P,sum.constraint){
  beta=exp(gamma)
  beta[gamma==0]=0
  M = Meila(F_Vec2Sym(beta))
  lambda = SetLambda(P, M,sum.constraint)
  D=length(gamma)
  #gradient with log transformation
  return((- F_Sym2Vec(P) + beta*(F_Sym2Vec(M) + lambda)))
}
F_NegLikelihood_Trans <- function(gamma, log.psi, P,sum.constraint){
  #gamma=gamma-mean(gamma)
  M = Meila(F_Vec2Sym(exp(gamma)))

  lambda = SetLambda(P, M,sum.constraint)

  suppressWarnings(
    res<-(-sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) )+
      log(SumTree(F_Vec2Sym(exp(gamma))))+
      lambda*(sum(exp(gamma))-sum.constraint/2))
  # cat("like val is... ",res," !!\n Detail: a=",
  # -sum(F_Sym2Vec(P) * (log(exp(gamma))+ F_Sym2Vec(log.psi))) ,"/b=",
  # log(SumTree(F_Vec2Sym(exp(gamma)))),"/c=",
  #   lambda*(sum(exp(gamma))-sum.constraint/2),"\n"
  # )
  #cat(paste0("\nSumTree=",SumTree(F_Vec2Sym(exp(gamma)))))


  return( res)
}
sum.constraint.inf<-function(p,min.order=300){
  round( p*(p-1)*10^(-min.order/(p-1)))+1
}
#########################################################################
SetLambda <- function(P, M,sum.constraint=1, eps = 1e-6, start=1){
  # F.x has to be increasing. The target value is 0

  F.x <- function(x){
    if(x!=0){
      sum.constraint - sum(P / (x+M))
    }else{
      sum.constraint - (2*sum(P[upper.tri(P)] / M[upper.tri(M)]))
    }
  }
  i=1
  x.min = ifelse(F.x(0) >0,-sort(F_Sym2Vec(M))[1]+1e-10,max(-1e-4, -min(F_Sym2Vec(M))/2));
  while(F.x(x.min)>0 && i<length(unique(M))){
    i=i+1
    x.min=-sort(F_Sym2Vec(M))[i]+1e-10
  }
  if(F.x(x.min)>0) stop("Could not set lambda.")
  x.max = start
  while(F.x(x.max)<0){x.max = x.max * 10}
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
  alpha.grid = (1:n)/n
  cond = Inf; a = 0
  while(cond > cond.tol & a<n){
    a = a+1
    psi.vec = -alpha.grid[a]*n*log(1 - F_Sym2Vec(CorY^2))/2;
    #psi.vec = psi.vec - mean(psi.vec)
    psi = F_Vec2Sym(exp(psi.vec))
    lambda = svd(psi)$d
    cond = min(abs(lambda))/max(abs(lambda))

  }
  alpha = alpha.grid[a]

  psi.vec = -alpha.grid[a]*n*log(1 - F_Sym2Vec(CorY^2))/2;
  #  if(sum(is.na(psi.vec))!=0) browser()
  psi.vec = psi.vec - mean(psi.vec)
  psi = F_Vec2Sym(exp(psi.vec))
  return(list(psi=psi, alpha=alpha))
}

#########################################################################
#' Update the beta weights according to the gradient ascent
#'
#' @param beta.init Initial beta weight matrix
#' @param psi Psi matrix, filled with ratios of bivariate probabilities over marginals, which can in the Gaussian
#' wase be deduced from the correlation matrix.
#' @param maxIter Maximum number of iterations
#' @param eps  Precision parameter controlling the convergence of weights beta
#' @param unlinked An optional vector of nodes which are not linked with each other
#' @param sum.weights Sum constraint for the weight matrix
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
#' @importFrom stats optim
#' @examples n=30
#' p=10
#' S=5
#' Y=data_from_scratch("tree",p=p,n=n)$data
#' beta = matrix(1/10,10,10)
#' psi=Psi_alpha(cor(Y), n)$psi
#' FitEM = FitBeta(beta.init=beta, psi=psi, maxIter = 6, sum.weights=sum(beta))
FitBeta <- function(beta.init, psi, maxIter=50, eps = 1e-6, unlinked=NULL,sum.weights){
  beta.tol = 1e-4
  beta.min = 1e-16
  beta.old = beta.init
  log.psi = log(psi+(psi==0))
  iter = 0
  p=nrow(beta.init)
  logpY = rep(0, maxIter)
  beta.diff = diff.loglik = 2 * eps

  stop=FALSE
  # min.val=max((-(p-2)*log(p)+round(log(.Machine$double.xmin))+10)/(p-1), -10)
  # max.val=min((-(p-2)*log(p)+round(log(.Machine$double.xmax))-10)/(p-1), 10)
  min.val= (-(p-2)*log(p)+round(log(.Machine$double.xmin))+10)/(p-1)
  max.val=(-(p-2)*log(p)+round(log(.Machine$double.xmax))-10)/(p-1)

  while (  ((beta.diff > eps) || (diff.loglik>eps) ) && iter < maxIter && !stop){
    iter = iter+1
    P=Kirshner(W=beta.old*psi)
    init=F_Sym2Vec(beta.old)

    #gradient ascent in log scale
    gamma_init=log(init+(init==0))
    gamma_init[init==0]=0

    gamma = stats::optim(gamma_init, F_NegLikelihood_Trans, gr=F_NegGradient_Trans,method='L-BFGS-B',
                         log.psi, P,sum.weights, control=list(trace=0,maxit=500,  pgtol=1e-2, factr=1e+8),
                         lower=rep(min.val, p*(p-1)/2),upper=rep(max.val, p*(p-1)/2))$par

    # cat(round(mean(gamma),2))
    beta=exp(gamma)
    beta[which(beta< beta.min)] = beta.min # numerical zeros
    beta=F_Vec2Sym(beta)
    if(!is.null(unlinked)) beta[unlinked, unlinked]=0
    logpY[iter] = -F_NegLikelihood(F_Sym2Vec(beta),log.psi,P,sum.weights)
    #if(iter==2)browser()
    beta.diff = max(abs(beta.old-beta))
    #cat(paste0("/ ",round(beta.diff,6),"\n"))
    beta.old = beta
    diff.loglik=logpY[iter]-logpY[iter-1]
    if(iter > 1){diff.loglik =  abs(diff.loglik)}else{diff.loglik=1}
  }
  logpY = logpY[1:iter]
  return(list(edges_weight=beta, logpY=logpY,maxIter=iter))
}

#########################################################################
#' Core computing function
#'
#' @param PLN.Cor Either an object resulting from the use of the `PLN` function from package `PLNmodels`, or an estimation of Gaussian data correlation matrix.
#' @param n Number of samples, required if a correlation matrix is supplied
#' @param maxIter Maximum number of iterations for EMtree
#' @param unlinked An optional vector of nodes which are not linked with each other
#' @param random.init A boolean for trying a random initialization of the EM
#' @param cond.tol Tolerance parameter for the conditioning of psi matrix
#' @param eps  Precision parameter controlling the convergence of weights beta
#' @param verbatim Talks if set to TRUE
#' @param plot Plots likelihood if set to TRUE
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
#' EMtree(PLN.Cor=PLN_Y,verbatim=TRUE, plot=TRUE)


EMtree<-function(PLN.Cor, n=NULL,  maxIter=30, unlinked=NULL , random.init=FALSE,
                 cond.tol=1e-10, eps = 1e-3,
                 verbatim=TRUE, plot=FALSE){
  T1<-Sys.time()
  if(inherits(PLN.Cor, "PLNfit")){
    CorY=cov2cor(PLN.Cor$model_par$Sigma)
    n=PLN.Cor$n
  }else if(inherits(PLN.Cor, "matrix") & nrow(PLN.Cor) == ncol(PLN.Cor)){
    CorY = PLN.Cor
  }else{
    stop("PLN.Cor must be a PLN object or a squarred gaussian correlation matrix")
  }
  p=ncol(CorY)


  # beta.init with adaptative mean value depending on the network dimensions
  # sum.weights=sum.constraint.inf(p)
  # mean.val=sum.weights/(p*(p-1))
  mean.val=exp((-(p-2)*log(p))/(p-1))
  beta.init = matrix(mean.val, p, p);  diag(beta.init)=0
  if(!is.null(unlinked)) beta.init[unlinked, unlinked]=0 #unliked nodes

  if(random.init){ # try different starting points for this EM
    beta.init = matrix(runif(n=p*p, min=0.9*mean.val,max=1.1*mean.val ), p,p)
    beta.init=t(beta.init)%*%beta.init/2
    diag(beta.init)=0
  }
  sum.weights=round(sum(beta.init))+1
  alpha.psi = Psi_alpha(CorY, n, cond.tol=cond.tol)
  psi = alpha.psi$psi

  lambda = svd(Laplacian(beta.init*psi)[-1,-1])$d
  cond = min(abs(lambda))/max(abs(lambda))
  adapt=FALSE
  if(cond<1e-12 || min(lambda)<=0){
    adapt=TRUE
    if(verbatim)cat("Adapting conditioning tolerance... ")
  }
  #  browser()
  while(cond<1e-12 || min(lambda)<=0){ # set the tempering parameter alpha, for a good conditioning of the psi matrix
    cond.tol=10*cond.tol
    alpha.psi = Psi_alpha(CorY, n, cond.tol=cond.tol)
    psi = alpha.psi$psi
    lambda = svd(Laplacian(beta.init*psi)[-1,-1])$d
    cond = min(abs(lambda))/max(abs(lambda))
  }
  if(adapt && verbatim)  cat(paste0("Final value:",cond.tol))
  FitEM = FitBeta(beta.init=beta.init, psi=psi, maxIter = maxIter,eps = eps,
                  unlinked=unlinked, sum.weights=sum.weights)

  beta=FitEM$edges_weight
  P = Kirshner(beta*psi)
  time<-difftime(Sys.time(),T1)

  if(plot){
    g<-tibble(p=FitEM$logpY) %>% rowid_to_column() %>%
      ggplot(aes(rowid,p))+geom_point()+geom_line()+theme_minimal()+labs(x="Iter",y="Likelihood")
    print(g)
  }

  if(verbatim){
    cat("\nConvergence took",round(time,2), attr(time, "units")," and ", FitEM$maxIter," iterations.")
  }
  return(c(list(edges_prob=P, norm.cst = SumTree(beta), timeEM=time),FitEM))

}

#################################################################################
#' Resampling procedure for  edges probability
#'
#' @param counts Data of observed counts with dimensions n x p, either a matrix, data.frame or tibble.
#' @param covar_matrix matrix of covariates, should have the same number of rows as the count matrix.
#' @param unlinked An optional vector of nodes which are not linked with each other
#' @param O Matrix of offsets, with dimension n x p
#' @param v The proportion of observed data to be taken in each sub-sample. It is the ratio (sub-sample size)/n
#' @param S Total number of wanted sub-samples.
#' @param maxIter Maximum number of EMtree iterations at each sub-sampling.
#' @param cond.tol Tolerance for the psi matrix.
#' @param eps Precision parameter controlling the convergence of weights beta
#' @param cores Number of cores, can be greater than 1 if data involves less than about 32 species.
#' @param init boolean: should the resmapling be carried out with different initial points (TRUE), or with different initial data (FALSE)
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
ResampleEMtree <- function(counts,covar_matrix=NULL  , unlinked=NULL, O=NULL,
                           v=0.8, S=1e2, maxIter=30, cond.tol=1e-10,eps=1e-3,cores=3, init=FALSE){
  cat("Computing ",S,"probability matrices with", cores, "core(s)... ")
  t1=Sys.time()
  counts=as.matrix(counts)
  n = nrow(counts);  p = ncol(counts)
  P = p * (p - 1) / 2 ; V = round(v * n)
  Pmat = matrix(0, S, P)
  #- offsets and covariates
  if(is.null(O)){ O=matrix(1, n, p)}
  if(is.null(covar_matrix)){#default intercept
    X=matrix(1,nrow=n,ncol=1)
  }else{X=as.matrix(covar_matrix)}
  #- parallel computation of S fits of new_EMtree
  suppressWarnings(
    PLNfit <- PLNmodels::PLN(counts ~ -1  + offset(log(O)) + ., data=data.frame(X), control = list("trace"=0))
  )
  obj<-parallel::mclapply(1:S,function(b){
    set.seed(b)
    if(init){
      inf<-EMtree( PLNfit,unlinked,n=n, maxIter=maxIter, cond.tol=cond.tol,verbatim=TRUE,eps=eps,
                   plot=FALSE, random.init = TRUE)[c("edges_prob","maxIter","timeEM")]
    }else{
      sample = sample(1:n, V, replace = F)

      #inception
      M.sample=PLNfit$var_par$M[sample,]
      S2.sample=PLNfit$var_par$S2[sample,]
      CorY=cov2cor(t(M.sample)%*%M.sample+diag(colSums(S2.sample)))
      try({
        inf<-EMtree( CorY,unlinked,n=n, maxIter=maxIter, cond.tol=cond.tol,verbatim=TRUE,eps=eps,
                     plot=FALSE, random.init=TRUE)[c("edges_prob","maxIter","timeEM")]
      }, silent=TRUE)
      if(!exists("inf")) inf=NA #depending on the sample drawn, it is possible that computation fail
      # because of bad conditioning of the Laplacian matrix of the weights beta.
      # This can happen especially when using the "unlinked" parameter.
    }
    return(inf)
  }, mc.cores=cores)
  bad_samples=which(do.call(rbind, lapply(obj, length))!=3)
  time=difftime(Sys.time(), t1)
  cat(round(time,2),  attr(time, "units"),"\n")
  if(length(bad_samples)!=0){
    cat(length(bad_samples), " failed samples.\n")
    obj=obj[-bad_samples]
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
#' @param unlinked An optional vector of nodes which are not linked with each other
#' @param Pt Probability threshold for every sub-sample
#' @param v The proportion of observed data to be taken in each sub-sample. It is the ratio (sub-sample size)/n
#' @param S Number of desired sub-samples
#' @param maxIter Maximum number of iterations of EMtree
#' @param cond.tol Tolerance for the psi matrix.
#' @param cores Number of cores, can be greater than 1 if data involves less than about 32 species.
#'
#' @return a tibble in a suitable format for graphical purposes. It has four columns describing each edge of each model: the endpoints ("node1", "node2"),
#'  the model it corresponds to ("model") and its weight ("value").
#' @export
#' @importFrom parallel mclapply
#' @importFrom purrr map
#' @importFrom tidyr gather  unnest
#' @importFrom dplyr mutate filter
#' @importFrom tibble tibble rownames_to_column
#' @importFrom  stats formula model.matrix
#' @examples
#'n=30
#'p=10
#'S=3
#'Y = data_from_scratch(type="tree",p=p,n=n)$data
#'X = data.frame(rnorm(n),rbinom(n,1,0.7))
#'ComparEMtree(Y,X,models=list(1,2,c(1,2)),m_names=list("1","2","both"),Pt=0.3,S=S, cores=1)
ComparEMtree <- function(counts, covar_matrix, models, m_names, O=NULL,unlinked=NULL, Pt=0.1, v=0.8,
                         S=1e2, maxIter=50, cond.tol=1e-14,cores=3){

  Stab.sel<- lapply(seq_along(models),function(x){
    cat("model",m_names[[x]],": ")
    covariates=colnames(covar_matrix)[models[[x]]]
    chaine=paste0("~",paste(covariates,collapse="+")) #includes the intercept
    formule=stats::formula(chaine)
    matcovar=stats::model.matrix(formule,covar_matrix)
    Pmat=ResampleEMtree(counts,matcovar,unlinked= unlinked,O=O, v=v, S=S, maxIter, cond.tol=cond.tol,cores=cores)$Pmat
    cat("\n")
    return(Pmat)
  })
  p=ncol(counts)

  mat<-data.frame(freq_selec(Stab.sel[[1]],Pt)) # the first element is initialized
  allNets<-tibble(P = list(mat), mods =m_names )  %>%
    mutate(P=purrr::map( seq_along(P), function(x) {
      df<-freq_selec(Stab.sel[[x]],Pt=Pt)
      df[lower.tri(df, diag = TRUE)]<-NA
      df<-data.frame(df)
      colnames(df)<-1:ncol(df)
      df
    })) %>%
    mutate(P = purrr::map(P,~rownames_to_column(.) %>%
                            gather(key, value , -rowname) %>%filter(!is.na(value))
    ),
    mods=unlist(mods)
    ) %>%
    tidyr::unnest(cols = c(P))
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
  E=ncol(Pmat) #number of edges
  p=sqrt(2*E+(1/4))-1/2
  if(is.null(Pt)) Pt= 2/p + 0.1
  return(F_Vec2Sym(colMeans(1*(Pmat>Pt))))
}


