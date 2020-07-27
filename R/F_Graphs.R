

#' Plots a natwork
#'
#' @param adj_matrix graph adjacency matrix
#' @param title graph title
#' @param size size of nodes
#' @param curv edges curvature
#' @param width maximum width for the edges
#' @param alpha if TRUE, sets to transparent the edges non-linked to nodes with high betweenness
#' @param filter_deg selects nodes with a higher degree than filter_deg
#' @param nb sets the number of nodes selected by thresholding the beetweenness scores
#' @param layout optional ggraph layout.
#' @param nodes_label optional labels for nodes.
#' @param pal optional palette.
#' @param seed optional seed for graph reproductibility.
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
#' @examples adj_matrix= SimCluster(p=30,k=3,dens=0.4, r=50)
#' draw_network(adj_matrix,"Cluster graph", layout="fr",curv=0.1)
draw_network<-function(adj_matrix,title="", size=4, curv=0,width=1, alpha=FALSE, filter_deg=FALSE,nb=1,
                       layout=NULL,nodes_label=NULL,pal="#31374f",  groupes=NULL){
  p=nrow(adj_matrix)
  if(is.null(nb)) nb=round(p/8,0)
  binary=FALSE
  if(is.null(nodes_label)){ nodes_label=1:p ; bool=TRUE}else{bool=FALSE}

  if(sum(unique(adj_matrix))==1) binary=TRUE

  res<- as_tbl_graph(adj_matrix, directed=FALSE) %>% activate(edges) %>%
    mutate(btw.weights=ifelse(weight==0,0,log(1+1/weight))) %>%
    activate(nodes) %>%
    mutate(  btw=centrality_betweenness(weights=btw.weights),
            bool_btw=(btw>sort(btw, decreasing = TRUE)[nb]),bool_deg=(centrality_degree()>0),
            deg=centrality_degree(), title=title, name=nodes_label
    )
  if(!is.null(groupes)){
    res<-res %>% mutate(groupes=as.factor(groupes))
  }
  if(bool){
    res<-res %>% mutate(label=ifelse(bool_btw,name,""))
  }else{
    res<-res %>% mutate(label=ifelse(bool_deg,name,""))
  }
  if(filter_deg) res <- res %>% activate(nodes) %>% filter(deg!=0)
  res<-res %>%
    activate(edges)  %>%
    filter(weight !=0) %>%
    mutate(neibs=edge_is_incident(which(.N()$bool_btw)), title=title)


  pal_edges <-  ifelse(is.null(pal), viridisLite::viridis(5, option = "C")[c(3,2,4,1)], pal)
 if(!is.null(groupes)){
   pal_nodes<-c("#31374f","#adc9e0","#e7bd42")
   res<-res %>%
     activate(nodes)  %>%
     mutate(finalcolor=groupes)
 }else{
   pal_nodes<-c("black","goldenrod1")
   res<-res %>%
     activate(nodes)  %>%
     mutate(finalcolor=bool_btw)
 }
  set_graph_style(family="sans")
  layout = ifelse(is.null(layout), "circle", layout)
  g=res %>%
    ggraph(layout = layout)
  if(alpha){
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
    scale_size_manual(values = c(2,5,6))+
    geom_node_text(aes(label = label), color = "black", size = size) +#,nudge_x = 0.3
    labs(title = title) + theme(plot.title = element_text(hjust = 0.5))

  if(!binary){ g=g+
    scale_edge_width_continuous(range=c(0.1,width))
  }else{
    g=g+
      scale_edge_width_continuous(range=c(0,width))
  }

  return(list(G=g,graph_data=res))

}




#' Plots networks from several models
#'
#' @param allNets tibble resulting of ComparEMtree()
#' @param curv edges curvature
#' @param width maximum width for the edges
#' @param alpha if TRUE, sets to transparent the edges non-linked to nodes with high betweenness
#' @param Ft Frequency threshold
#' @param seed optional seed for graph reproductibility
#' @param nb sets the number of nodes selected by thresholding the beetweenness scores
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
#' @examples
#'n=30
#'p=10
#'S=3
#'Y=data_from_scratch("tree",p=p,n=n)$data
#'X = data.frame(rnorm(n),rbinom(n,1,0.7))
#'threemodels=ComparEMtree(Y,X,models=list(1,2,c(1,2)),
#'m_names=list("1","2","both"),Pt=0.3,S=S, cores=1)
#'compar_graphs(threemodels)
compar_graphs<-function(allNets, curv=0.2, width=1, alpha=TRUE,Ft=0,
                        nodes_label=NULL,seed=123, nb=3, layout="circle", base_model=NULL){

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
        mutate( importance=centrality_degree(),btw=centrality_betweenness(weights = btw.weights),boolbtw=(btw>sort(btw, decreasing = TRUE)[nb]),
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

  set.seed(seed)
  plot<- spliT$P %>%
    reduce(bind_graphs) %>%
    activate(nodes) %>%
    mutate(mod=factor(mod,levels=mods),x=rep(lay$x,nbmod),y=rep(lay$y,nbmod)) %>%
    activate(edges) %>% filter(weight>Ft) %>%
    ggraph(layout="nicely")
  if(alpha){
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

