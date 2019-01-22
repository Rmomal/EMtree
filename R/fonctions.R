#  DICTIONNAIRE


####################################################
#g??n??ration de graph
# @graphER donne une matrice sym??trique al??atoire de taille n*n obtenue par
#          bernoulli de param p
####################################################
makeSymm <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

graphER<-function(n,p){
  matAdj<-matrix(rbinom(n*n, 1, p), ncol = n, nrow = n)
  matAdj<-makeSymm(matAdj)
  diag(matAdj)<-c(rep(0,n))
  return(matAdj)
}
#require(igraph)
# library(igraph)
# igraph_opt<-function(x, default = NULL){
#   if (missing(default)) {
#     get_config(paste0("igraph::", x), .igraph.pars[[x]])
#   }
#   else {
#     get_config(paste0("igraph::", x), default)
#   }
# }
# sample_gnp<-function(n, p, directed = FALSE, loops = FALSE){
#   type <- "gnp"
#   type1 <- switch(type, gnp = 0, gnm = 1)
#   on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
#   res <- .Call("R_igraph_erdos_renyi_game", as.numeric(n),
#                as.numeric(type1), as.numeric(p), as.logical(directed),
#                as.logical(loops), PACKAGE = "igraph")
#   if (igraph_opt("add.params")) {
#     res$name <- sprintf("Erdos renyi (%s) graph", type)
#     res$type <- type
#     res$loops <- loops
#     res$p <- p
#   }
#   res
# }
# SpannTree<-function(d){
#   g <- sample_gnp(d, 1)
#   g_mst <- mst(g)
#   adj<-as.matrix(as_adj(g_mst))
#   return(adj)
# }
library(vegan)
SpannTree <- function(p){
  W = matrix(runif(p^2), p, p); W = W + t(W)
  Tree = spantree(W)
  G = matrix(0, p, p)
  invisible(sapply(1:length(Tree$kid), function(i){G[i+1, Tree$kid[i]] <<- 1}))
  G = G + t(G)
  return(G)
}
erdos<-function(d,prob){
  G = matrix(rbinom(d^2, 1, prob), d, d)
  G[lower.tri(G, diag=T)] = 0; G = G+t(G)
}
SimCluster <- function(p, k, d, r){
  # k = nb clusters, d = graph density, r = within/between ratio connection probability
  beta = d / (r / k + (k - 1) / k)
  alpha = r * beta
  while (alpha > 1) {
    r = .9 * r
    beta = d / (r / k + (k - 1) / k)
    alpha = r * beta
  }
  Z = t(rmultinom(p, 1, rep(1 / k, k)))
  ZZ = Z %*% t(Z)
  diag(ZZ) = 0
  ZZ = F_Sym2Vec(ZZ)
  G = F_Vec2Sym(rbinom(p * (p - 1) / 2, 1, alpha * ZZ + beta * (1 - ZZ)))
  #gplot(G, gmode='graph', label=1:p, vertex.col=Z%*%(1:k),
  #     main=paste0(round(r,2),"// alpha ",round(alpha,2),"// beta ",round(beta,2)))
  
  return(G)
}
# d=10
# par(mfrow=c(2,2))
# SimCluster(d,3,5/d,2)
# SimCluster(d,3,5/d,5)
# SimCluster(d,3,5/d,9)
# SimCluster(d,3,5/d,15)
generator_graph<-function(d = 20, graph = "tree", g = NULL, prob = NULL, dens=0.3, vis = FALSE,
                          verbose = TRUE,r=5){
  gcinfo(FALSE)
  if (verbose)
    #cat("Generating data from the multivariate normal distribution with the",
    #     graph, "graph structure....")
    if (is.null(g)) {
      g = 1
      if (graph == "hub" || graph == "cluster") {
        if (d > 40)
          g = ceiling(d/20)
        if (d <= 40)
          g = 2
      }
    }
  theta = matrix(0, d, d)
  if (graph == "cluster") {
    #browser()
    
    theta<-SimCluster(d,3,dens,r) #prob=5/d ?
    
  }
  if (graph == "scale-free") {
    
    out = .C("SFGen", dd0 = as.integer(2), dd = as.integer(d),
             G = as.integer(theta), PACKAGE = "huge")
    theta = matrix(as.numeric(out$G), d, d)
  }
  if(graph=="tree"){
    theta<-SpannTree(d)
  }
  if(graph=="erdos"){
    theta<- erdos(d=d,p=prob)
    if(sum(theta)<4){
      while(sum(theta)<4){
        theta<- erdos(d=d,p=prob)
      }
    }
  }
  
  #browser()
  if (verbose)
    cat("done.\n")
  rm(vis, verbose)
  gc()
  return(theta = Matrix(theta, sparse = TRUE))
}
generator_param<-function(G,v=1){
  cste = 1
  omega = diag(rep(cste, ncol(G))) + G*v
  while (min(eigen(omega)$values) < 1e-6){
    cste = 1.1*cste
    omega = diag(rep(cste, ncol(G))) + G*v
  }
  #browser()
  sigma = cov2cor(solve(omega))
  #omega = solve(sigma)
  sim=list(sigma=sigma,omega=omega,cste=cste)
  return(sim)
}
####################################################
# construction de matrice de pr??cision
# @curseur utilise le r??sultat de g??n??ration de matrice al??atoire symm??trique par
#         une fonction donn??e et permet de faire varier la diagonale par g
# @findG permet d'obtenir la plus petite valeure de g possible pour avoir une
#         matrice d??finie positive
# @Omega utilise le g optimal pour g??n??rer la matrice de pr??cision voulue
####################################################
# librairies :library(matrixcalc)


precision<-function(g,fonction,argms){
  G<-fonction(argms[1], argms[2])
  Omega<-G+ g * diag(rowSums(G))
  return(Omega)
}

findG2<-function(seq,fonction,argms){
  liste_min_vp<-unlist(lapply(seq,function(x) min(eigen(precision(x,fonction,argms))$values)))
  res<-seq[which(liste_min_vp>0)[1]]
  if (is.na(res)){
    stop("Inappropriate sequence")
  }else{
    return(res)
  }
}

Omega<-function(fonction,argms,seq){
  g<-findG2(seq,fonction,argms)
  Omega<-precision(g,fonction,argms)
  return(Omega)
}

simMatrix<-function(n){
  #abs(rnorm(n*n))
  Sigma<-matrix(runif(n*n), ncol = n, nrow = n)
  Sigma<-makeSymm(Sigma)
  return(Sigma)
}

# exemple : Omega(graphER,c(20,0.3),seq(2,5,by=0.2),20)

####################################################
# Simulation et comparaison du graph avec l'originel
# @simu permet de simuler des variables ind??pendantes selon une loi normale
#       multi-vari??e de mat de pr??cision omega
# @roc utilise le package ROCR et permet de comparer les matrices omega et omega_hat. Cette
#       fonction est utile pour un plot direct ou superposition de plot faciles.
# @fun.auc.ggplot calcule UNE courbe ROC ?? partir de deux matrices. Donne les statistiques
#       de sens, spec et auc. Joli ggplot.
####################################################

# library(ROCR)
# library(mixtools)
# library(glasso)
# library(RColorBrewer)


estim_reg<-function(Y,G){
  d<-ncol(Y)
  n<-nrow(Y)
  if(n>d){ 
    S = P  = matrix(0, d,d)
    for (j in 1:d){
      #  browser()
      LM = lm(Y[, j] ~ -1 + Y[, -j])
      P[j, -j] = summary(LM)$coef[, 4]
      S[j, -j] = summary(LM)$coef[, 3]
    }
    # hist(as.vector(P), breaks=p)
    # boxplot(as.vector(S) ~ as.vector(Gdiag))
    # qqnorm(as.vector(S[which(Gdiag==0)]))
    res<-c(2*sum(P>.5), sum(G==0))
  }else{
    res<-c("not enough obs",sum(G==0))
  }
  return(res)
}


# simu<-function(omega,n.vectors){
#   sigma<-makeSymm(solve(omega))
#   mu<-c(rep(0,ncol(sigma)))
#   sim<-rmvnorm(n.vectors,mu,sigma)
#   return(sim)
# }

roc<-function(mat_rho,omega){
  indices_nuls<-which(omega==0)
  labels<-matrix(1,nrow=nrow(omega),ncol=ncol(omega))
  labels[indices_nuls]<-0
  pred<-prediction(as.vector(mat_rho),as.vector(labels))
  perf <- performance(pred,"tpr","fpr")
  invisible(perf)
}
fun.auc.ggplot <- function(pred, obs, title,points){
  nvar<-ncol(obs)
  obs[which(abs(obs)<1e-16)]<-0
  indices_nuls<-which(obs==0)
  label<-matrix(1,nrow=nrow(obs),ncol=ncol(obs))
  label[indices_nuls]<-0
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),as.vector(label[upper.tri(label)]))
  obs<-as.vector(label[upper.tri(label)])
  # Run the AUC calculations
  ROC_perf <- performance(prediction,"tpr","fpr")
  ROC_sens <- performance(prediction,"sens","spec")
  ROC_auc <- performance(prediction,"auc")
  
  # Make plot data
  plotdat <- data.frame(FP=ROC_perf@x.values[[1]],TP=ROC_perf@y.values[[1]],CUT=ROC_perf@alpha.values[[1]],POINT=NA)
  
  if(max(points)<max(plotdat$CUT[is.finite(plotdat$CUT)])){
    points<-points[-which(points>max(plotdat$CUT[is.finite(plotdat$CUT)]))]
  }
  
  #plotdat[unlist(lapply(points,function(x){which.min(abs(plotdat$CUT-x))})),"POINT"]<- points
  
  # Plot the curve
  rows<-which(!is.na(plotdat$POINT))
  while(length(rows)>20){
    rows<-rows[seq(1,length(rows),2)]
  }
  
  ggplot(plotdat, aes(x=FP,y=TP)) +
    geom_abline(intercept=0,slope=1,linetype="dashed",color="grey") +
    geom_line() +
    geom_point(data=plotdat[rows,], aes(fill=POINT),color="royalblue3", pch=20, size=3) +
    geom_text(data=plotdat[rows,],  label=round(plotdat$POINT[rows],5), hjust=-0.2, vjust=1.5) +
    scale_fill_gradientn("Threhsold Cutoff",colours=brewer.pal(5,"RdPu"),guide=FALSE) +
    scale_x_continuous("False Positive Rate", limits=c(0,1)) +
    scale_y_continuous("True Positive Rate", limits=c(0,1)) +
    geom_polygon(aes(x=X,y=Y), data=data.frame(X=c(0.77,1,1,0.77),Y=c(0,0,0.33,0.33)), fill="white") +
    annotate("text",x=0.97,y=0.30,label=paste0("Nnods = ",nvar),hjust=1) +
    annotate("text",x=0.97,y=0.25,label=paste0("Nedges = ",sum(obs==1)),hjust=1) +
    annotate("text",x=0.97,y=0.20,label=paste0("Nvoids = ",sum(obs==0)),hjust=1) +
    annotate("text",x=0.97,y=0.15,label=paste0("AUC = ",round(ROC_auc@y.values[[1]],digits=2)),hjust=1) +
    annotate("text",x=0.97,y=0.10,label=paste0("Sens = ",round(mean(as.data.frame(ROC_sens@y.values)[,1]),digits=2)),hjust=1) +
    annotate("text",x=0.97,y=0.05,label=paste0("Spec = ",round(mean(as.data.frame(ROC_sens@x.values)[,1]),digits=2)),hjust=1) +
    theme(legend.position="none", plot.title=element_text(vjust=2)) +
    labs(title=title)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
}

diagnostic.auc.sens.spe <- function(pred, obs,stat="auc",method){
  nvar<-ncol(obs)
  obs[which(abs(obs)<1e-16)]<-0
  
  indices_nuls<-which(obs==0)
  label<-matrix(1,nrow=nrow(obs),ncol=ncol(obs))
  label[indices_nuls]<-0
  prediction<-prediction(as.vector(pred[upper.tri(pred)]),
                         as.vector(label[upper.tri(label)]))
  obs<-as.vector(label[upper.tri(label)])
  # Run the AUC calculations
  ROC_sens <- performance(prediction,"sens","spec")
  ROC_auc <- performance(prediction,"auc")
  ROC_precision<-performance(prediction,"prec")
  ROC_recal<-performance(prediction,"rec")
  precrec<-data.frame(ROC_precision@y.values,ROC_recal@y.values,method)[-1,]
  colnames(precrec)<-c("prec","rec","method")
  #plot(ROC_precision)
  res<-switch(stat,"auc"=round(ROC_auc@y.values[[1]],digits=3),
              "sens"=round(mean(as.data.frame(ROC_sens@y.values)[,1]),digits=3),
              "spec"=round(mean(as.data.frame(ROC_sens@x.values)[,1]),digits=3))
  return(list(res,precrec))
}

build_crit<-function(path, nbgraph,x,variable,method, crit="auc"){ #method="_spiecResid" 
  if(variable=="n"){
    obs<-readRDS(paste0(path,"/Sets_param/Graph",nbgraph,".rds"))$omega
  }else{
    obs<-readRDS(paste0(path,"/Sets_param/Graph",nbgraph,"_",x,".rds"))$omega  
  } 
  # if(method=="_treeggm_" || method=="_oracle"|| method=="_one_step_" || method=="_mint_"){
  #   pred<- readRDS(paste0(path,"/Scores/Graph",nbgraph,method,x,".rds"))[[colonne]]
  # }else{
  #   pred<- readRDS(paste0(path,"/Scores/Graph",nbgraph,method,x,".rds"))}
  
    pred<- readRDS(paste0(path,"/Scores/",method,"_graph",nbgraph,"_",x,".rds"))
#mehtod: EMmarg, EMCond, OneMarg, OneCond
    res<-switch(crit,"auc"=diagnost_auc(obs,pred), "precrec"=diagnost_precrec(obs,pred),
                "precrecPool"=vec_obs_pred(obs,pred))
  return(res)
}

vec_obs_pred<-function(obs, pred){
  nvar<-ncol(obs)
  obs[which(abs(obs)<1e-16)]<-0
  indices_nuls<-which(obs==0)
  label<-matrix(1,nrow=nrow(obs),ncol=ncol(obs))
  label[indices_nuls]<-0
  
  vec_pred<-as.vector(pred[upper.tri(pred)])
  vec_obs<-as.vector(label[upper.tri(label)])
  
  return(list(vec_pred,vec_obs))
}
diagnost_precrec<-function(obs, pred){
  obs_pred<-vec_obs_pred(obs,pred)
  prediction<-prediction(obs_pred[[1]],obs_pred[[2]])
  ROC_precision<-performance(prediction,"prec")
  ROC_recal<-performance(prediction,"rec")
  return(list(ROC_precision@y.values,ROC_recal@y.values))
}
diagnost_auc<-function(obs, pred){
  obs_pred<-vec_obs_pred(obs,pred)
  prediction<-prediction(obs_pred[[1]],obs_pred[[2]])
  # Run the AUC calculations
  ROC_auc <- performance(prediction,"auc")
  res<-round(ROC_auc@y.values[[1]],digits=3)
  return(res)
}
####################################################
# Exploration des valeurs critiques de rho
# @tableau3D cr??er le tableau ?? trois dimension compos?? des matrices omega_hat pour
#           des rho diff??rents
# @monotonie permet l'??tude de la monotonie du ph??nom??ne d'extinction d'une ar??te (case
#           de la matrice omega_hat). Renvoie les indices en trois dimensions des incoh??rences.
# @mat_rho donne, selon la valeur du param??tre minmax, la valeur min ou la valeur max
#         de rho qui annule une case de omega_hat. Retourne la heatmap non r??ordonn??e
#         des valeurs critiques, ainsi que leur densit??, en ??chelle log.
####################################################
# install.packages(c("grid","gridGraphics","gridExtra"))
# library(grid)
# library(gridGraphics)
# library(gridExtra)

tableau3D<-function(samples,seq_rho){
  n<-dim(var(samples))[1]
  tab <- array( dim=c(n,n,length(seq_rho)))
  #nb_nuls<-data.frame(matrix(ncol=2,nrow=length(seq_rho)))
  for(rho in seq_rho){
    
    omega_hat<-glasso(var(samples),rho=rho,penalize.diagonal = FALSE)$wi
    tab[,,which(seq_rho==rho)]<-omega_hat
    #nb_nuls[which(seq_rho==rho),"nuls"]<-length(which(omega_hat==0))*100/length(omega_hat)
    #nb_nuls[which(seq_rho==rho),"rho"]<-rho
  }
  # g<-ggplot(nb_nuls,aes(log10(rho),nuls))+
  #   geom_line()+
  #   labs(y="% of nul elements")
  # print(g)
  return(tab)
}

defaut_monotonie<-function(tab,seq_rho){
  non_monotone<-0
  L<-length(seq_rho)-1
  infos_non_mono<-data.frame(matrix(ncol=3))
  colnames(infos_non_mono)<-c("ligne","colonne","rho")
  for(ligne in 1:n){
    for(colonne in ligne:n){
      for(rho in 1:L){
        if(tab[ligne,colonne,rho]==0 & tab[ligne,colonne,rho+1]!=0){
          non_monotone<-non_monotone+1
          infos_non_mono[non_monotone,1]<-ligne
          infos_non_mono[non_monotone,2]<-colonne
          infos_non_mono[non_monotone,3]<-seq_rho[rho]
        }
      }
    }
  }
  return(infos_non_mono)
}

mat_rho<-function(tab,seq_rho,minmax){
  results<-NA*tab[,,1]
  results[which(tab[,,1]==0)]<-0
  sequence<-switch(minmax,"min"=1:length(seq_rho),"max"=(length(seq_rho)-1):1)
  #NB : dans l'autre sens on cherche le dernier nul
  for(index.rho in sequence){
    #browser()
    indices_nuls<-switch(minmax,"min"=which(tab[,,index.rho]==0),
                         "max"=which(tab[,,index.rho]!=0))
    if(length(indices_nuls)!=0){
      indices_remplis<-which(!is.na(results))
      if(length(indices_remplis)!=0){
        indices_a_remplir<-setdiff(indices_nuls,indices_remplis) #indices nuls - les indices remplis
        
      }else{
        indices_a_remplir<-indices_nuls
      }
      results[indices_a_remplir]<-switch(minmax,"min"=seq_rho[index.rho],"max"=seq_rho[index.rho+1])
    }
  }
  diag(results)<-NA
  # par(mfrow=c(1,2))
  # heatmap.2(results,symm=TRUE,dendrogram="none",main=minmax,
  #           Colv=FALSE,Rowv=FALSE,density.info="none",keysize=1.5,trace="none",col=brewer.pal(9,"Blues")[3:9])
  # # q<-grab_grob()
  # g<-ggplot(data.frame(Penalty=as.vector(results[upper.tri(results)])),aes(log10(Penalty)))+
  #   geom_density()+
  #   labs(title=minmax,hjust=0.5)
  # grid.arrange(q,g, ncol=2, clip=TRUE)
  return(results)
}


#s = cov(data)
####################################################
# fonctions d'inf??rence
####################################################
get.lambda.l1 <- function(S) {
  r.max <- max(0,max(abs(S-diag(diag(S)))))
  return(r.max)
}
inf_glasso_MB<-function(X){
  S<-cov(X)
  d<-ncol(X)
  log.lambda.min <- -5
  log.lambda.max <- log(get.lambda.l1(S))
  log.lambda <- seq(log.lambda.min, log.lambda.max, length.out = d*2)
  MB.res <- lapply(exp(log.lambda),
                   function(lambda)
                     glasso(S, lambda, trace = FALSE, approx = FALSE, penalize.diagonal = FALSE))
  adjmat.array <- simplify2array(Map("*",
                                     exp(log.lambda),
                                     lapply(MB.res, function(x){ (abs(x$wi)>0)*1})
  )
  )
  # Let us replace each edge by the  largest Glasso lambda where it disappears (or a sum related to this)
  K.score <- apply(adjmat.array,c(1,2),sum)
  K.score <- K.score / max(K.score)
  return(K.score)
}
#install_github("zdk123/SpiecEasi")
#library(SpiecEasi)
inf_spieceasi<-function(Y){
  inf<-spiec.easi(Y, icov.select = FALSE, nlambda = 50, verbose = FALSE)
  
  #browser()
  adjmat.array <- simplify2array(Map("*",
                                     inf$lambda,
                                     lapply(inf$est$path, function(x) {
                                       (abs(x) > 0) * 1
                                     })))
  K.score <- Reduce("+",adjmat.array)
  K.score <- K.score / max(K.score)
}
inf_spiecresid<-function(Y,covar){
  U<-clr(Y)
  model<-lm(U~covar)
  inf<-inf_spieceasi(model$residuals)
  return(inf)
}

clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

clr.matrix <- function(x.f, mar=2, ...) {
  apply(x.f, mar, clr, ...)
}



####################################################
# Consid??rations gaphiques
# @paramNet cr??er des param??tres graphiques par d??faut pour les graphs
####################################################
#librairaies :library(igraph, RColorBrewer)

paramNet<-function(net,pal,label){
  deg <- degree(net, mode="out")
  V(net)$size <- deg*2+5
  V(net)$color=pal[3]
  E(net)$color=pal[8]
  V(net)$name <-label
  E(net)$curved=.1
  E(net)$width<-rep(2, ecount(net))
  V(net)$frame.color=pal[3]
  return(net)
}
matrice_adj<-function(precision,seuil){
  adj<-unlist(apply(precision,1, function(x) ifelse(abs(x)<seuil,0,1)))
  return(adj)
}

net_from_matrix<-function(precision,seuil, boucles){
  matrice<-matrice_adj(precision,seuil)
  net <- graph_from_adjacency_matrix(matrice, mode="upper", diag = boucles)
  return(net)
}
net_from_weight<-function(precision,seuil, boucles){
  net <- graph_from_adjacency_matrix(precision, mode="upper", diag = boucles,weighted=TRUE)
  return(net)
}
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  # return gtable invisibly
  invisible(combined)
}

