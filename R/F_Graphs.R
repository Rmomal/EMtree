

#' Plots a network
#'
#' @param adj_matrix graph adjacency matrix
#' @param title graph title
#' @param label_size size of labels
#' @param curv edges curvature
#' @param width maximum width for the edges
#' @param shade if TRUE, shades the edges unlinked to nodes with high betweenness
#' @param remove_isolated Removes isolated nodes if TRUE
#' @param btw_rank betweenness rank +1 of highlighted nodes. If set to 1, none are highlighted.
#' @param layout optional layout from `ggraph` or `ìgraph` package.
#' @param stored_layout optional data.frame of point positions from a previous graph, possibly created using ggraph::create_layout()
#' @param nodes_label optional labels for nodes.
#' @param nodes_size size of nodes, possibility to specify a size per group
#' @param pal_edges optional palette for edges
#' @param pal_nodes optional palette for nodes
#' @param node_groups optional vector separating the nodes into groups
#' @param edge_groups optional matrix specifying groups of edges
#' @param legend optional boolean for generating a legend when groups are provided
#' @param leg.pos optional legend position ("bottom", "top", "right" or "left")
#' @param title.family font family of the title
#' @param title.size font size of the title
#' @return \itemize{
#' \item{G} {the network as a ggplot2 object, with highlighted high betweenness nodes}
#' \item{graph_data}{ data needed for plotting the network}
#' \item{legend}{ ggplot2 information about the legend}
#' }
#' @export
#' @importFrom tidygraph .N as_tbl_graph activate centrality_betweenness centrality_degree edge_is_incident
#' @importFrom ggraph ggraph set_graph_style geom_edge_arc scale_edge_alpha_manual scale_edge_width_continuous geom_node_text geom_node_point  scale_edge_colour_manual create_layout
#' @importFrom tibble tibble rownames_to_column
#' @importFrom ggplot2 aes theme labs scale_color_manual scale_size_manual ggplot_gtable ggplot_build element_text
#' @importFrom dplyr mutate filter
#' @importFrom tibble as_tibble
#' @importFrom gridExtra grid.arrange
#' @examples
#' set.seed(1)
#' adj_matrix= SimCluster(p=20,k=2,dens=0.5, r=10)
#' groups=sample(3,20, replace=TRUE)
#' draw_network(adj_matrix,"Cluster graph", layout="stress",
#' shade=TRUE, btw_rank=3)$G
#' draw_network(adj_matrix,"Cluster with groups of nodes",
#' layout="stress", shade=TRUE, btw_rank=3,
#' node_groups=groups, legend=TRUE,title.family="mono",title.size=12)$G
draw_network<-function(adj_matrix,title="", label_size=4, curv=0,width=1,
                       shade=FALSE, remove_isolated=FALSE,btw_rank=2,
                       layout=NULL,stored_layout=NULL,nodes_label=NULL,
                       nodes_size=c(2,5),pal_edges=NULL, pal_nodes=NULL,
                       node_groups=NULL, edge_groups=NULL,legend=FALSE,
                       leg.pos="bottom", title.family="sans",title.size=12){
  adj_matrix<-as.matrix(adj_matrix)
  p<-nrow(adj_matrix) ; binary<-FALSE

  if(is.null(nodes_label)){ nodes_label<-1:p ; nonames<-TRUE
  }else{nonames<-FALSE}
  if(sum(unique(adj_matrix))==1) binary<-TRUE
 if(!is.null(edge_groups)){
   edge_groups<-ToVec(edge_groups)[ToVec(adj_matrix!=0)]
 }else{
   edge_groups<-ToVec(matrix(1, p,p))[ToVec(adj_matrix!=0)]
 }
  edge_groups<-as.factor(edge_groups)
  # edges width
  min.width<-ifelse(binary,0,0.1)
  #edges colour
if(is.null(pal_edges)){
  ncol_e<-length(unique(edge_groups))
  if(ncol_e<5){
    pal_edges<-c("#31374f","#1F968BFF","#B8DE29FF","de7065ff")[1:ncol_e]
  }else{stop("Please fill the pal_edges parameter")}

}
  #betweenness computations

  res<- as_tbl_graph(adj_matrix, directed=FALSE) %>% activate(edges) %>%
    mutate(btw.weights=ifelse(weight==0,0,log(1+1/weight))) %>%
    activate(nodes) %>%
    mutate( btw=centrality_betweenness(weights=btw.weights),
            bool_btw=(btw>sort(btw, decreasing = TRUE)[btw_rank]),
            bool_deg=(centrality_degree()>0),
            deg=centrality_degree(), title=title, name=nodes_label )
  if(!is.null(node_groups)) res<-res %>%
    mutate(node_groups=as.factor(node_groups))

  if(nonames){
    res<-res %>% mutate(label=ifelse(bool_btw,name,""))
    #if no names supplied, only print important labels
  }else{
    res<-res %>% mutate(label=name)
  }
  if(remove_isolated) res <- res %>% activate(nodes) %>% filter(deg!=0)
  res<-res %>%
    activate(edges)  %>%
    filter(weight !=0) %>%
    mutate(neibs=edge_is_incident(which(.N()$bool_btw)), title=title)

  # define nodes color

  sensitive_nodes<-unlist(res %>%  activate(nodes)  %>%
                           dplyr::select(bool_btw) %>%
                            tibble::as_tibble())
  if(!is.null(node_groups)){
    if(is.null(pal_nodes)) pal_nodes<-c("#31374f","#adc9e0","#e7bd42")
    if(sum(sensitive_nodes)==0) sensitive_nodes<-node_groups
  }else{
    node_groups<-sensitive_nodes
    if(is.null(pal_nodes)) pal_nodes<-c("#31374f","#e7bd42")
  }

  #draw graph
  set_graph_style(title.family)
  layout <- ifelse(is.null(layout), "nicely", layout)

  if(is.null(stored_layout)){
    finallayout<-create_layout(res,layout=layout)
    g<-res %>%
      ggraph(graph = finallayout)
  }else{
    finallayout<-stored_layout[,1:2]
    g<-res %>%
      ggraph(layout = finallayout)
  }


  if(shade){ #shading
    g<-g+
      geom_edge_arc(aes(edge_width=weight, alpha=neibs,color = edge_groups),
                    strength = curv, show.legend = FALSE) +
      scale_edge_alpha_manual(values=c(0.2,1))
  }else{ g<-g+
    geom_edge_arc(aes(edge_width=weight,color = edge_groups), strength = curv,
                  show.legend = FALSE)
  }

  g<-g+
    geom_node_point(aes(color = node_groups, size = as.factor(sensitive_nodes)),
                    show.legend =c(colour=legend,size=FALSE)) +
    scale_edge_colour_manual(values = pal_edges) +
    scale_color_manual("",values = pal_nodes)+
    scale_size_manual(values = nodes_size)+
    geom_node_text(aes(label = label), color = "black", size = label_size) +
    labs(title = title) + theme(plot.title = ggplot2::element_text(
      hjust=0.5,
      family = title.family,
      size = title.size ),
      legend.position=leg.pos)+
    scale_edge_width_continuous(range=c(min.width,width))
leg<-NULL
  if(!is.null(node_groups) & legend ){#add legend if groups provided
    tmp<-ggplot(data.frame(node_groups=as.factor(node_groups),
                           row1=adj_matrix[1,], row2=adj_matrix[2,]),
               aes(row1, row2, color=node_groups))+
      geom_point()+scale_color_manual("",values=pal_nodes)+
      theme(legend.position=leg.pos)
    tmp <- ggplot_gtable(ggplot_build(tmp))
    leg <- tmp$grobs[[which(lapply(tmp$grobs, function(x) x$name) == "guide-box")]]

  }

  return(list(G=g,graph_data=res,legend=leg))
}

utils::globalVariables(c("edges", "weight", "nodes", "btw.weights", "btw",
                         "bool_btw","name","deg","neibs","label","row1","row2"))

