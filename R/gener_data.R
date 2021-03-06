

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
  W <- matrix(runif(p^2), p, p); W <- W + t(W)
  Tree <- spantree(W)
  G <- matrix(0, p, p)
  invisible(
    lapply(seq_along(Tree$kid), function(i){G[i+1, Tree$kid[i]] <<- 1})
    )
  G <- G + t(G)
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
  G <- matrix(rbinom(p^2, 1, prob), p,p)
  G[lower.tri(G, diag=TRUE)] <- 0; G <- G+t(G)
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
#' @importFrom stats rmultinom
#' @examples SimCluster(10,2,0.5, 10)
SimCluster <- function(p, k, dens, r){
  beta <- dens / (r / k + (k - 1) / k)
  alpha <- r * beta
  while (alpha > 1) {
    r <- .9 * r
    beta <- dens / (r / k + (k - 1) / k)
    alpha <- r * beta
  }
  Z <- t(stats::rmultinom(p, 1, rep(1 / k, k)))
  groupe<-Z%*%1:k
  Z <- Z %*% t(Z)
  diag(Z) <- 0
  ZZ <- ToVec(Z)
  G <- ToSym(rbinom(p * (p - 1) / 2, 1, alpha * ZZ + beta * (1 - ZZ)))
  res<-G
  return(res)
}
#' Simulate several types of graphs
#'
#' @param p number of nodes
#' @param graph type of graph, among "tree","scale-free","cluster" and "erdos"
#' @param dens graph density (for cluster graphs) or edges probability (for Erdös-Renyi graphs)
#' @param r within/between ratio connection probability (needed for cluster graphs)
#' @param k number of groups in cluster graphs
#'
#' @return the adjacency matrix, in sparse format
#' @export
#' @importFrom  Matrix Matrix
#' @importFrom huge huge.generator
#' @examples generator_graph(p=10,graph="tree")
generator_graph<-function(p = 20, graph = "tree", dens=0.3, r=2, k=3){
  theta <- matrix(0, p, p)
  if (graph == "cluster") {
    theta<-SimCluster(p,k=k,dens,r)
  }
  if (graph == "scale-free") {
    theta <- huge::huge.generator(d=p,graph="scale-free",verbose = FALSE)$theta
  }
  if(graph=="tree"){
    theta<-SpannTree(p)
  }
  if(graph=="erdos"){
    theta<- erdos(p=p,prob=dens)
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
#'  \item{lambda: }{the constant that was needed to load the diagonal of omega and ensure its positive definiteness}
#' }
#' @export
#' @importFrom stats rbinom
#' @examples G<-generator_graph(p=10,graph="tree")
#' generator_param(G=G)
#'

generator_param<-function(G,signed=FALSE,v=0.01){
  lambda <- 1
  p<-ncol(G)
  sumlignes<-rowSums(matrix(G,p,p))
  if(sum(sumlignes==0)!=0) sumlignes[sumlignes==0]<-0.1
  D<-diag(sumlignes+v)
  if(signed){
    Gsign <- ToSym(ToVec(G * matrix(2*stats::rbinom(p^2, 1, .3)-1, p, p)))
    omega <- lambda*D + Gsign
    while(min(eigen(omega)$values) < 1e-10 & lambda<1e3){
      lambda <- 1.1*lambda
      omega <- lambda*D + Gsign
    }
  }else{
    omega <- lambda*D + G
    while (min(eigen(omega)$values) < 1e-10){
      lambda <- 1.1*lambda
      omega <-lambda*D + G
    }
  }
  sigma <- solve(omega)
  sim<-list(sigma=sigma,omega=omega,lambda=lambda)
  return(sim)
}
#' Simulate count data under the Poisson log-Normal model
#'
#' @param Sigma Covariance matrix of the normal hidden layer of parameters
#' @param covariates a data.frame or matrix containing data covariates. If not NULL, defines the
#' number of  simulated rows.
#' @param n number of rows to simulate
#' @param norm should the parameters be normalized ?
#'
#' @return  \itemize{
#' \item{if norm=FALSE}{ Y: the simulated counts}
#' \item{if norm=TRUE}{ \itemize{
#'       \item{Y:}{simulated counts}
#'       \item{U:}{ the normalized Gaussian parameters}}}}
#' @export
#' @importFrom  mvtnorm rmvnorm
#' @importFrom stats as.formula model.matrix cov2cor rpois runif
#' @examples G<-generator_graph(p=10,graph="tree")
#' sigma<-generator_param(G=G)$sigma
#' generator_PLN(as.matrix(sigma))
generator_PLN<-function(Sigma,covariates=NULL, n=50, norm=FALSE){
  p<-ncol(Sigma)
  if(!is.null(covariates)){
    n<-nrow(covariates)
    c<-ncol(covariates)

    string<-paste0("~", paste(colnames(covariates), collapse=" + "))
    formula<-stats::as.formula(string)
    m<- stats::model.matrix(formula,covariates)

    mc<-ncol(m)-1
    beta<-matrix(stats::runif(p*mc),mc,p)
    prod<-(matrix(m[,-1],n , mc) %*% beta) +2
  }else{
    prod<-2 # constant for signal
  }
  if(norm){
    D<-diag(Sigma)
    R<-stats::cov2cor(Sigma)
    U<- mvtnorm::rmvnorm(n, rep(0,p), R)
    matsig<-(matrix(rep(sqrt(D),n),n,p, byrow = TRUE))
    Y <- matrix(stats::rpois(n*p, exp(U*matsig+prod )), n, p)
    sim<-list(Y=Y, U=U)
  }else{
    Z<- mvtnorm::rmvnorm(n, rep(0,nrow(Sigma)), Sigma)
    Y <- matrix(stats::rpois(n*p, exp(Z+prod )), n, p)
    sim<-Y
  }
  return(sim)
}
#data_from_scratch:

#' Simulates data under the PLN model with control on the dependency structure
#'
#' @param type type of graph, either "tree", "erdos", "cluster" or "scale-free"
#' @param p wanted number of columns (species)
#' @param n wanted number of rows (samples/observations)
#' @param r within/between connectivity ratio for cluster graphs
#' @param covariates a data.frame or matrix containing data covariates.
#' @param dens density of edges for cluster graphs or edge probability for Erdös graphs
#' @param k number of groups for the cluster structure.
#' @param norm should the Gaussian parameters be normalized ?
#' @param signed boolean: should the graph be composed of positive and negative partial correlations ?
#' @param v noise parameter of the precision matrix
#' @param draw boolean, plots the graph if set to TRUE
#' @return a list containing
#' \itemize{
#'  \item{data: }{simulated counts}
#'  \item{omega: }{the precision matrix}
#' }
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom tidygraph as_tbl_graph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point theme_graph
#' @examples set.seed(1)
#' p<-30
#' Y1<-data_from_scratch("tree",p=p,draw=TRUE)
#' str(Y1)
#' Y2<-data_from_scratch("cluster",p=p,r=10, dens=10/p, k=3,draw=TRUE)
data_from_scratch<-function(type, p=20,n=50, r=5, covariates=NULL,
                            dens=log(p)/p,k=3, norm=FALSE, signed=FALSE,
                            v=0.01,draw=FALSE){
  # make graph
  graph<- generator_graph(graph=type,p=p,dens=dens,r=r, k=k)
  param<-generator_param(G=as.matrix(graph),signed=signed,v=v)
  data<-generator_PLN(param$sigma,covariates,n, norm=norm)
  if(draw){
    g<-as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="stress")+
      geom_edge_link(color="gray70")+#"#31374f"
      geom_node_point(size=2, color="steelblue")+theme_graph()
    print(g)
  }
  return(list(data=data,omega= param$omega))
}
