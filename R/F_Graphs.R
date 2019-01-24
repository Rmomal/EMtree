

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
draw_network<-function(df,adj_covar, pal=NULL){
  res<- as_tbl_graph(df, directed=FALSE) %>%
    activate(nodes) %>%
    mutate( keyplayer = node_is_keyplayer(k=3), model=adj_covar, name=1:p,
            label=ifelse(keyplayer,name,"")) %>%
    activate(edges)  %>%
    filter(weight !=0) %>%
    mutate(neibs=edge_is_incident(which(.N()$keyplayer)), model=adj_covar)

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
compar_graphs<-function(allNets, alpha=TRUE){
  nbmod<-length(unique(allNets$models))

  spliT<-data.frame(allNets) %>%
    split(allNets$models) %>%
    tibble(P=map(.,function(x){
      model<-x$models[1]
      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(nodes) %>%
        mutate( importance=centrality_degree(),
                keyplayer = node_is_keyplayer(k=3), model=model) %>%
        activate(edges) %>% filter(value!=0) %>%
        mutate(neibs=edge_is_incident(which(.N()$keyplayer)), model=model) %>%
        activate(nodes) %>%
        mutate(label=ifelse(keyplayer,name,"")) #%>%
      #  filter(importance!=0)
    }))
  mods=unique(allNets$models)
  pal_edges <- viridisLite::viridis(5, option = "C")
  pal_nodes<-c("gray15","goldenrod1")
  lay<-create_layout(spliT$P[2][[1]],layout="circle")
  set_graph_style()

  plot<- spliT$P[1][[1]] %>%
    bind_graphs(spliT$P[2][[1]] )%>%
    bind_graphs(spliT$P[3][[1]] )%>%
    bind_graphs(spliT$P[4][[1]] )%>%
    activate(nodes) %>%
    mutate(model=factor(model,levels=mods),x=rep(lay$x,4),y=rep(lay$y,4)) %>%
    ggraph(layout="auto")
  if(alpha){
    plot<-plot+
      geom_edge_arc(aes(color=model, alpha=neibs),curvature=0.3,show.legend=FALSE)+
      scale_edge_alpha_manual(values=c(0.2,1))

  }else{plot<-plot+
    geom_edge_arc(aes(color=model),curvature=0.3,show.legend=FALSE)
  }
  plot+
    geom_node_point(aes(color=keyplayer, size=keyplayer), show.legend=FALSE)+
    scale_edge_colour_manual(values=pal_edges[c(4,2,1,3)], labels=mods)+
    scale_color_manual(values=pal_nodes)+
    scale_size_manual(values=c(1.5,6))+
    geom_node_text(aes(label = label),color="black")+
    facet_nodes(~model, scales="free",ncol=4)+
    th_foreground(border=FALSE)+
    theme(strip.background = element_rect(fill="white",color="white"),
          strip.text = element_text(color="black",size=12))

}
