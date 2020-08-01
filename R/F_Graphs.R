

#' Plots a network
#'
#' @param adj_matrix graph adjacency matrix
#' @param title graph title
#' @param size size of nodes
#' @param curv edges curvature
#' @param width maximum width for the edges
#' @param shade if TRUE, shades the edges unlinked to nodes with high betweenness
#' @param filter_deg selects nodes with a higher degree than filter_deg
#' @param btw_rank betweenness rank +1 of highlighted nodes. If set to 1, none are highlighted.
#' @param layout optional ggraph layout.
#' @param nodes_label optional labels for nodes.
#' @param nodes_size size of nodes, possibility to specify a size per group
#' @param pal_edges optional palette for edges
#' @param pal_nodes optional palette for nodes
#' @param groupes optional vector seperating the nodes into groupes
#'
#' @return \itemize{
#' \item{G} {the network as a ggplot2 object, with highlighted high betweenness nodes}
#' \item{graph_data}{ data needed for plotting the network}
#' }
#' @export
#' @importFrom tidygraph .N as_tbl_graph activate centrality_betweenness centrality_degree edge_is_incident
#' @importFrom ggraph ggraph set_graph_style geom_edge_arc scale_edge_alpha_manual scale_edge_width_continuous geom_node_text geom_node_point  scale_edge_colour_manual
#' @importFrom tibble tibble rownames_to_column
#' @importFrom ggplot2 aes theme labs scale_color_manual scale_size_manual
#' @importFrom dplyr mutate filter
#' @importFrom viridisLite viridis
#' @importFrom tibble as_tibble
#' @examples adj_matrix= SimCluster(20,2,0.4, 10)
#' draw_network(adj_matrix,"Cluster graph", layout="fr", shade=TRUE)
draw_network<-function(adj_matrix,title="", size=4, curv=0,width=1, shade=FALSE, filter_deg=FALSE,btw_rank=2,
                       layout=NULL,nodes_label=NULL,nodes_size=c(2,5),pal_edges=NULL, pal_nodes=NULL, groupes=NULL){
  p=nrow(adj_matrix) ; binary=FALSE
  if(is.null(nodes_label)){ nodes_label=1:p ; nonames=TRUE}else{nonames=FALSE}
  if(sum(unique(adj_matrix))==1) binary=TRUE
  # edges width
  min.width=ifelse(binary,0,0.1)
  #edges colour
  pal_edges <-  ifelse(is.null(pal_edges), viridisLite::viridis(5, option = "C")[c(3,2,4,1)], pal_edges)
  #betweenness computations
  res<- as_tbl_graph(adj_matrix, directed=FALSE) %>% activate(edges) %>%
    mutate(btw.weights=ifelse(weight==0,0,log(1+1/weight))) %>%
    activate(nodes) %>%
    mutate( btw=centrality_betweenness(weights=btw.weights),
            bool_btw=(btw>sort(btw, decreasing = TRUE)[btw_rank]),
            bool_deg=(centrality_degree()>0),
            deg=centrality_degree(), title=title, name=nodes_label )
  if(!is.null(groupes)) res<-res %>% mutate(groupes=as.factor(groupes))
  if(nonames){
    res<-res %>% mutate(label=ifelse(bool_btw,name,""))
  }else{
    res<-res %>% mutate(label=ifelse(bool_deg,name,""))
  }
  if(filter_deg) res <- res %>% activate(nodes) %>% filter(deg!=0)
  res<-res %>%
    activate(edges)  %>%
    filter(weight !=0) %>%
    mutate(neibs=edge_is_incident(which(.N()$bool_btw)), title=title)

# define nodes color
  if(!is.null(groupes)){
   if(is.null(pal_nodes)) pal_nodes<-c("#31374f","#adc9e0","#e7bd42")
 }else{
   groupes=unlist(res %>%
     activate(nodes)  %>% dplyr::select(bool_btw) %>% tibble::as_tibble())
   if(is.null(pal_nodes)) pal_nodes<-c("#31374f","#e7bd42")
 }
  res<-res %>%
    activate(nodes)  %>%
    mutate(finalcolor=groupes)

  #draw graph
  set_graph_style(family="sans")
  layout = ifelse(is.null(layout), "circle", layout)
  g=res %>% ggraph(layout = layout)
  if(shade){ #shading
    g<-g+
      geom_edge_arc(aes(edge_width=weight, alpha=neibs,color = title), strength = curv, show.legend = FALSE) +
      scale_edge_alpha_manual(values=c(0.2,1))
  }else{ g<-g+
    geom_edge_arc(aes(edge_width=weight,color = title), strength = curv, show.legend = FALSE)
  }
  g<-g+
    geom_node_point(aes(color = finalcolor, size = groupes), show.legend = FALSE) +
    scale_edge_colour_manual(values = pal_edges) +
    scale_color_manual(values = pal_nodes)+
    scale_size_manual(values = nodes_size)+
    geom_node_text(aes(label = label), color = "black", size = size) +#,nudge_x = 0.3
    labs(title = title) + theme(plot.title = element_text(hjust = 0.5))+
    scale_edge_width_continuous(range=c(min.width,width))

  return(list(G=g,graph_data=res))
}




#' Plots networks from several models
#'
#' @param allNets tibble resulting of ComparEMtree()
#' @param curv edges curvature
#' @param width maximum width for the edges
#' @param shade if TRUE, shades the edges non-linked to nodes with high betweenness
#' @param Ft Frequency threshold
#' @param btw_rank betweenness rank +1 of highlighted nodes. If set to 1, none are highlighted.
#' @param layout choice of layout for all networks among available layouts in `ggraph`. Default is "circle"
#' @param base_model choice of referent model for the layout construction
#' @param nodes_label optional label for the nodes
#'
#' @return \itemize{
#' \item{G}{ the collection of networks as a ggplot2 object, with highlighted high betweenness nodes}
#' \item{graph_data }{list of data needed to plot the networks}
#' }
#' @export
#' @importFrom tidygraph .N bind_graphs as_tbl_graph activate centrality_betweenness centrality_degree edge_is_incident
#' @importFrom ggraph  ggraph create_layout set_graph_style facet_nodes th_foreground geom_edge_arc scale_edge_alpha_manual scale_edge_width_continuous geom_node_text geom_node_point  scale_edge_colour_manual
#' @importFrom tibble tibble rownames_to_column
#' @importFrom ggplot2 element_rect element_text aes theme labs scale_color_manual scale_size_manual
#' @importFrom dplyr mutate filter
#' @importFrom purrr map reduce
#' @importFrom viridisLite viridis
#' @examples
#'n=30
#'p=10
#'S=3
#'Y=data_from_scratch("tree",p=p,n=n)$data
#'X = data.frame(rnorm(n),rbinom(n,1,0.7))
#'threemodels=ComparEMtree(Y,X,models=list(1,2,c(1,2)),
#'m_names=list("1","2","both"),Pt=0.3,S=S, cores=1)
#'compare_graphs(threemodels)
compare_graphs<-function(allNets, curv=0.2, width=1, shade=TRUE,Ft=0, nodes_label=NULL,
                         btw_rank=2, layout="circle", base_model=NULL){

   mods=unique(allNets$model)
  if(is.null(base_model)) base_model = allNets$model[1]
  nbmod<-length(mods)
  binary=(sum(unique(allNets$weight))==1)

if(is.null(nodes_label)){
  nb_sp = max(as.numeric(allNets$node2))
  nodes_label=1:nb_sp
}
  spliT<-data.frame(allNets) %>%
    base::split(allNets$model) %>%
    tibble(P=map(.,function(x){
      mod<-x$model[1]

      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>% filter(weight!=0) %>%
        mutate(btw.weights=log(1+1/weight)) %>%
        activate(nodes) %>%
        mutate( importance=centrality_degree(),btw=centrality_betweenness(weights = btw.weights),boolbtw=(btw>sort(btw, decreasing = TRUE)[btw_rank]),
                 mod=mod, names=nodes_label) %>%
        activate(edges) %>%
        mutate(neibs=edge_is_incident(which(.N()$boolbtw)), mod=mod) %>%
        activate(nodes) %>%
        mutate(label=ifelse(boolbtw,names,""))

    }))


  pal_edges <- viridisLite::viridis(5, option = "C")
  pal_nodes<-c("gray15","goldenrod1")
  lay<-create_layout(spliT$P[[base_model]],layout=layout)
  set_graph_style(family="sans")
  plot<- spliT$P %>%
    reduce(bind_graphs) %>%
    activate(nodes) %>%
    mutate(mod=factor(mod,levels=mods),x=rep(lay$x,nbmod),y=rep(lay$y,nbmod)) %>%
    activate(edges) %>% filter(weight>Ft) %>%
    ggraph(layout="nicely")
  if(shade){
    plot<-plot+
      geom_edge_arc(aes(edge_width=weight,color=mod, alpha=neibs),strength=curv,show.legend=FALSE)+
      scale_edge_alpha_manual(values=c(0.2,1))

  }else{plot<-plot+
    geom_edge_arc(aes(edge_width=weight,color=mod),strength=curv,show.legend=FALSE)
  }

  g= plot+
    geom_node_point(aes(color=boolbtw, size=boolbtw), show.legend=FALSE)+
    scale_edge_colour_manual(values=pal_edges[c(1,3,2,4)], labels=mods)+
    scale_color_manual(values=pal_nodes)+
    scale_size_manual(values=c(1.5,6))+
    geom_node_text(aes(label = label),color="black", repel = FALSE)+
    facet_nodes(~mod, scales="free",ncol=nbmod)+
    th_foreground(border=FALSE)+
    theme(strip.background = element_rect(fill="white",color="white"),
          strip.text = element_text(color="black",size=14))
  if(!binary){ g=g+
    scale_edge_width_continuous(range=c(0.1,width))
  }else{
    g=g+
      scale_edge_width_continuous(range=c(0,width))
  }
  return(list(G=g, graph_data=spliT))
}

