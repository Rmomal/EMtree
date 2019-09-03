

#' Title
#'
#' @param adj_matrix graph adjacency matrix
#' @param title graph title
#' @param size size of nodes
#' @param curv edges curvature
#' @param filter_deg selects nodes with a higher degree than filter_deg
#' @param layout optional igraph layout.
#' @param nodes_label optional labels for nodes.
#' @param pal optional palette.
#' @param seed optional seed for graph reproductibility.
#'
#' @return plots a network
#' @export
#' @import tidygraph dplyr ggraph ggplot2
#'
#' @examples adj_matrix= SimCluster(10,2,0.5, 0.1)
#' draw_network(adj_matrix,"Cluster graph")
draw_network<-function(adj_matrix,title, size=4, curv=0.3, filter_deg=FALSE,layout=NULL,nodes_label=NULL,pal=NULL,
                       seed=200){
  p=nrow(adj_matrix)
  nb=round(p/6,0)
  if(is.null(nodes_label)){ nodes_label=1:p ; bool=TRUE}else{bool=FALSE}

  res<- as_tbl_graph(adj_matrix, directed=FALSE) %>%
    activate(nodes) %>%
    mutate( keyplayer = node_is_keyplayer(k=3), btw=centrality_betweenness(),
            bool_btw=(btw>sort(btw, decreasing = TRUE)[nb]),bool_deg=(centrality_degree()>0),
            deg=centrality_degree(), title=title, name=nodes_label
    )
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
    pal_nodes<-c("gray15","goldenrod1")

    set_graph_style(family="sans")
    layout = ifelse(is.null(layout), "circle", layout)
    res %>%
      ggraph(layout = layout) +
      geom_edge_arc(aes(color = title), curvature = curv, show.legend = FALSE) +
      geom_node_point(aes(color = bool_btw, size = bool_btw), show.legend = FALSE) +
      scale_edge_colour_manual(values = pal_edges) +
      scale_color_manual(values = pal_nodes) +
      scale_size_manual(values = c(1.5, 6)) +
      geom_node_text(aes(label = label), color = "black", size = size) +
      labs(title = title) + theme(plot.title = element_text(hjust = 0.5))

}




#' Title
#'
#' @param allNets object resulting of ComparEMtree()
#' @param alpha if TRUE, sets edges non-linked to nodes with high betweenness transparent
#' @param seed optional seed for graph reproductibility
#' @param nb sets the number of nodes selected by thresholding the beetweenness scores
#'
#' @return plots a collection of networks
#' @export
#' @import tidygraph dplyr ggraph ggplot2 influenceR
#'
#' @examples
compar_graphs<-function(allNets, alpha=TRUE,seed=123, nb=3, pos=1){

  nbmod<-length(unique(allNets$models))

  spliT<-data.frame(allNets) %>%
    split(allNets$models) %>%
    tibble(P=map(.,function(x){
      model<-x$models[1]

      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>% filter(value!=0) %>%
        activate(nodes) %>%
        mutate( importance=centrality_degree(),btw=centrality_betweenness(),boolbtw=(btw>sort(btw, decreasing = TRUE)[nb]),
                keyplayer = node_is_keyplayer(k=3), model=model) %>%
        activate(edges) %>%
        mutate(neibs=edge_is_incident(which(.N()$boolbtw)), model=model) %>%
        activate(nodes) %>%
        mutate(label=ifelse(boolbtw,name,""))

    }))

  mods=unique(allNets$models)
  pal_edges <- viridisLite::viridis(5, option = "C")
  pal_nodes<-c("gray15","goldenrod1")
  lay<-create_layout(spliT$P[[pos]],layout="circle")
  set_graph_style()
set.seed(seed)
  plot<- spliT$P %>%
    reduce(bind_graphs) %>%
    activate(nodes) %>%
    mutate(model=factor(model,levels=mods),x=rep(lay$x,nbmod),y=rep(lay$y,nbmod)) %>%
    ggraph(layout="auto")
  if(alpha){
    plot<-plot+
      geom_edge_arc(aes(color=model, alpha=neibs),curvature=0.3,show.legend=FALSE)+
      scale_edge_alpha_manual(values=c(0.2,1))

  }else{plot<-plot+
    geom_edge_arc(aes(color=model),curvature=0.3,show.legend=FALSE)
  }

  plot+
    geom_node_point(aes(color=boolbtw, size=boolbtw), show.legend=FALSE)+
    scale_edge_colour_manual(values=pal_edges[c(1,2,4,3)], labels=mods)+
    scale_color_manual(values=pal_nodes)+
    scale_size_manual(values=c(1.5,6))+
    geom_node_text(aes(label = label),color="black")+
    facet_nodes(~model, scales="free",ncol=nbmod)+
    th_foreground(border=FALSE)+
    theme(strip.background = element_rect(fill="white",color="white"),
          strip.text = element_text(color="black",size=14))

}

