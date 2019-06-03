

#' Title
#'
#' @param df
#' @param adj_covar
#' @param pal
#'
#' @return
#' @export
#'
#' @examples
draw_network<-function(df,adj_covar, pal=NULL,seed=200){
  set.seed(seed)
  res<- as_tbl_graph(df, directed=FALSE) %>%
    activate(nodes) %>%
    mutate( keyplayer = node_is_keyplayer(k=3), model=adj_covar, name=1:p,
            label=ifelse(keyplayer,name,"")) %>%
    activate(edges)  %>%
    filter(weight !=0) %>%
    mutate(neibs=edge_is_incident(which(.N()$keyplayer)), model=adj_covar)
#browser()
  pal_edges <-  ifelse(is.null(pal), viridisLite::viridis(5, option = "C")[c(3,2,4,1)], pal)
  pal_nodes<-c("gray15","goldenrod1")

  set_graph_style()
  res %>%
    ggraph(layout="circle")+
    geom_edge_arc(aes(color=adj_covar),curvature=0.3,show.legend=FALSE)+
    geom_node_point(aes(color=keyplayer, size=keyplayer), show.legend=FALSE)+
    scale_edge_colour_manual(values=pal_edges)+
    scale_color_manual(values=pal_nodes)+
    scale_size_manual(values=c(1.5,6))+
    geom_node_text(aes(label = label),color="black")+
    labs(title=adj_covar)+theme(plot.title = element_text(hjust = 0.5))
}




#' Title
#'
#' @param allNets
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
compar_graphs<-function(allNets, alpha=TRUE,seed=123, nb=3){

  nbmod<-length(unique(allNets$models))

  spliT<-data.frame(allNets) %>%
    split(allNets$models) %>%
    tibble(P=map(.,function(x){
      set.seed(seed)
      model<-x$models[1]
      #browser()
      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>% filter(value!=0) %>%
        activate(nodes) %>%
        mutate( importance=centrality_degree(),btw=centrality_betweenness(),boolbtw=(btw>sort(btw, decreasing = TRUE)[nb]),
                keyplayer = node_is_keyplayer(k=3), model=model) %>%
        activate(edges) %>%
        mutate(neibs=edge_is_incident(which(.N()$boolbtw)), model=model) %>%
        activate(nodes) %>%
        mutate(label=ifelse(boolbtw,name,"")) #%>%
      #  filter(importance!=0)
    }))
#  browser()
  mods=unique(allNets$models)
  pal_edges <- viridisLite::viridis(5, option = "C")
  pal_nodes<-c("gray15","goldenrod1")
  lay<-create_layout(spliT$P[4][[1]],layout="circle")
  set_graph_style()

  plot<- spliT$P[1][[1]] %>%
    bind_graphs(spliT$P[3][[1]] )%>%
    bind_graphs(spliT$P[2][[1]] )%>%
    bind_graphs(spliT$P[4][[1]] )%>%
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

