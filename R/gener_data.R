

##############
# DATA
##############
#data_from_scratch: generate data under the PLN model with a certain type of dependency structure,
#                   and draw the latter structure if asked.
# type: type of graph, either "tree", "erdos", "cluster" or "scale-free"
# p: wanted number of columns (species)
# covar: covariates
# prob: edge probability for erdos graphs
# dens: density of edges for cluster graphs
# r: within/between connectiviy ratio for cluster graphs

data_from_scratch<-function(type, p=20, r=5, covar,prob=log(p)/p,dens=log(p)/p, draw=FALSE){
  # make graph
  graph<- generator_graph(graph=type,d=p,prob=prob,dens=dens,r=r)
  param<-generator_param(as.matrix(graph))
  data<-generator_PLN(param$sigma,covar)[[1]]
  if(draw){ as_tbl_graph(as.matrix(graph)) %>%
      ggraph(layout="kk")+
      geom_edge_link()+
      geom_node_point(size=3, color="blue")
  }
  return(list(data=data,omega= param$omega))
}


