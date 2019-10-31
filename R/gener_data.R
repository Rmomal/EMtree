

##############
# DATA
##############
#' Simulates a spanning tree graph from `vegan`package
#'
#' @param p number of nodes
#'
#' @return the adjacency matrix
#' @export
#' @importFrom vegan spantree
#' @examples SpannTree(10)
SpannTree <- function(p){
  W = matrix(runif(p^2), p, p); W = W + t(W)
  Tree = spantree(W)
  G = matrix(0, p, p)
  invisible(sapply(1:length(Tree$kid), function(i){G[i+1, Tree$kid[i]] <<- 1}))
  G = G + t(G)
  return(G)
}
#' Simulates an Erdös-Renyi graph
#'
#' @param p number of nodes
#' @param prob probability for an edge, which is the same for all edges
#'
#' @return the adjacency matrix
#' @export
#'
#' @examples erdos(10,0.3)
erdos<-function(p,prob){
  G = matrix(rbinom(p^2, 1, prob), p,p)
  G[lower.tri(G, diag=T)] = 0; G = G+t(G)
  return(G)
}

#' Simulates a cluster graph
#'
#' @param p number of nodes
#' @param k number of groups in the cluster
#' @param dens graph density
#' @param r within/between ratio connection probability
#'
#' @return the adjacency matrix of size p x p
#' @export
#'
#' @examples SimCluster(10,2,0.5, 10)
SimCluster <- function(p, k, dens, r){
  beta = dens / (r / k + (k - 1) / k)
  alpha = r * beta
  while (alpha > 1) {
    r = .9 * r
    beta = dens / (r / k + (k - 1) / k)
    alpha = r * beta
  }
  Z = t(rmultinom(p, 1, rep(1 / k, k)))
  Z = Z %*% t(Z)
  diag(Z) = 0
  ZZ = F_Sym2Vec(Z)
  G = F_Vec2Sym(rbinom(p * (p - 1) / 2, 1, alpha * Z + beta * (1 - Z)))
  return(G)
}
#' Simulate several types of graphs
#'
#' @param p number of nodes
#' @param graph type of graph, among "tree","scale-free","cluster" and "erdos"
#' @param prob edges probability (for erdös-renyi graphs)
#' @param dens graph density (for cluster graphs)
#' @param r within/between ratio connection probability (needed for cluster graphs)
#'
#' @return the adjacency matrix, in sparse format
#' @export
#' @importFrom  Matrix Matrix
#' @importFrom huge huge.generator
#' @examples generator_graph(p=10,graph="tree")
generator_graph<-function(p = 20, graph = "tree", prob = 0.1, dens=0.3, r=5){
  theta = matrix(0, p, p)
  if (graph == "cluster") {
    theta<-SimCluster(p,3,dens,r)
  }
  if (graph == "scale-free") {
    theta = huge.generator(d=p,graph="scale-free")$theta

  }
  if(graph=="tree"){
    theta<-SpannTree(p)
  }
  if(graph=="erdos"){
    theta<- erdos(p=p,prob=prob)
  }
  return(theta = Matrix(theta, sparse = TRUE))
}
#' Generate simulation parameters from a graph adjacency matrix
#'
#' @param G adjacency matrix
#' @param signed boolean: should the graph be composed of positive and negative partial correlations ?
#' @param v parameter controlling the noise on the precision matrix
#'
#' @return a list containing
#' \itemize{
#'  \item{sigma: }{the covariance matrix}
#'  \item{omega: }{the precision matrix}
#'  \item{cste: }{the constant that was needed to load the diagonal of omega and esure its positive definitiveness}
#' }
#' @export
#'
#' @examples G=generator_graph(p=10,graph="tree")
#' generator_param(G=G)
#'

generator_param<-function(G,signed=FALSE,v=0.01){
  lambda = 1
  p=ncol(G)
  sumlignes=rowSums(matrix(G,p,p))
  D=diag(sumlignes+v)
  if(signed){
    Gsign = F_Vec2Sym(F_Sym2Vec(G * matrix(2*rbinom(p^2, 1, .3)-1, p, p)))
    omega = lambda*D + Gsign
    while(min(eigen(omega)$values) < 1e-10 & lambda<1e3){
      lambda = 1.1*lambda
      omega = lambda*D + Gsign
    }
    print(lambda)
  }else{
    print(dim(G))
    omega = lambda*D + G
    while (min(eigen(omega)$values) < 1e-10){
      lambda = 1.1*lambda
      omega =lambda*D + G
    }
  }
  sigma = cov2cor(solve(omega))
  sim=list(sigma=sigma,omega=omega,cste=lambda)
  return(sim)
}
#' Simulate count data under the Poisson log-Normal model
#'
#' @param Sigma Covariance matrix of the normal hidden layer of parameters
#' @param covariates a data.frame or matrix containing data covariates. If not NULL, defines the
#' number of  simulated rows.
#' @param n number of rows to simulate
#'
#' @return  Y: the simulated counts
#' @export
#' @importFrom  mvtnorm rmvnorm
#' @examples G=generator_graph(p=10,graph="tree")
#' sigma=generator_param(G=G)$sigma
#' generator_PLN(as.matrix(sigma))
generator_PLN<-function(Sigma,covariates=NULL, n=50){
  p<-ncol(Sigma)
  if(!is.null(covariates)){
    n<-nrow(covariates)
    c<-ncol(covariates)

    string<-paste0("~", paste(colnames(covariates), collapse=" + "))
    formula<-as.formula(string)
    m<- model.matrix(formula,covariates)[,-1]

    mc<-ncol(m)
    beta<-matrix(runif(p*mc),mc,p)
    prod=m %*% beta
  }else{
    prod=0
  }
  Z<- rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
  Y = matrix(rpois(n*p, exp(Z+prod )), n, p)
  return(Y)
}
#data_from_scratch:

#' generates data under the PLN model with a certain type of dependency structure, and draws the structure.
#'
#' @param type type of graph, either "tree", "erdos", "cluster" or "scale-free"
#' @param p wanted number of columns (species)
#' @param n wanted number of rows (samples/observations)
#' @param r within/between connectiviy ratio for cluster graphs
#' @param covariates a data.frame or matrix containing data covariates.
#' @param prob edge probability for erdos graphs
#' @param dens density of edges for cluster graphs
#' @param draw boolean, plots the graph if set to TRUE
#'
#' @return a list containing
#' \itemize{
#'  \item{data: }{simulated counts}
#'  \item{omega: }{the precision matrix}
#' }
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom tidygraph as_tbl_graph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#' @examples data_from_scratch("tree",p=10,draw=TRUE)
data_from_scratch<-function(type, p=20,n=50, r=5, covariates=NULL,prob=log(p)/p,dens=log(p)/p, draw=FALSE){
  # make graph
  graph<- generator_graph(graph=type,p=p,prob=prob,dens=dens,r=r)
  param<-generator_param(G=as.matrix(graph))
  data<-generator_PLN(param$sigma,covariates,n)
  if(draw){ as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="kk")+
      geom_edge_link()+
      geom_node_point(size=2, color="blue")
  }
  return(list(data=data,omega= param$omega))
}
