# PLNtree for Barents fish data

library(PLNmodels)
library(ade4)
library(tidyverse)
source('/Users/raphaellemomal/these/pack1/R/codes/FunctionsMatVec.R')
source('/Users/raphaellemomal/these/pack1/R/codes/FunctionsTree.R')
source('/Users/raphaellemomal/these/pack1/R/codes/FunctionsInference.R')

# Data
data.dir = '/Users/raphaellemomal/these/pack1/R/Data_RM/'
data(baran95)
data.name = 'Baran'
Data = list(
  count = as.matrix(baran95$fau),
  covariates = baran95$plan,
  species.names = baran95$species.names
)
Data$covariates$date = as.factor(Data$covariates$date)
Data$covariates$site = as.factor(Data$covariates$site)
colnames(Data$count) = baran95$species.names
n = nrow(Data$count)
p = ncol(Data$count)
species = colnames(Data$count)

# Algo parms
freq.sel.thres = 0.80
iter.max = 1e2
B = 5e2
VEM.fit = T
Tree.fit = T
Tree.res = T
v=0.8
# Fit VEM-PLN with 1 = no covariates, 2 = all covariates
if (VEM.fit) {
  #  VEM = list()
  
  #  VEM[[2]] = PLN(Data$count ~ Data$covariates$date)
  
  # VEM[[4]] = PLN(Data$count ~ Data$covariates$date + Data$covariates$site)
  # VEM[[5]] = PLN(Data$count ~ Data$covariates$date * Data$covariates$site)
  # save(VEM, file = paste0(data.dir, data.name, '-VEM.Rdata'))
  obj<- mclapply(1:B,function(b){
    #  cat('\n', b, '')
    set.seed(b);V = round(v*n)
    sample = sample(1:n, V, replace = F)
    Y.sample = Data$count[sample,]
    X.sample = Data$covariates[sample,]
    O.sample = O[sample,]

    VEM<-tibble(null=diag(PLN( Y.sample ~ 1)$model_par$Sigma), 
                site= diag(PLN( Y.sample ~  X.sample$site)$model_par$Sigma), graph=b)
    return(VEM)
  },mc.cores=1)
}
VEM<-do.call(rbind,obj)
VEM %>% mutate(diff=null-site) %>% group_by(graph) %>% summarise(med_diff=median(diff),sd_diff=sd(diff)) %>% 
  summary(med_diff)

load(paste0(data.dir, data.name, '-VEM.Rdata'))
M = length(VEM)

# Fit TreeGGM
if (Tree.fit) {
  EMtree = list()
  sapply(1:M, function(m) {
    EMtree[[m]] <<-
      TreeGGM(cov2cor(VEM[[m]]$model_par$Sigma),
              n,
              step = 'FALSE',
              maxIter = 500)
  })
  save(EMtree, file = paste0(data.dir, data.name, '-EMtree.Rdata'))
  
  sapply(1:M, function(m) {
    plot(
      F_Sym2Vec(EMtree[[m]]$P),
      F_Sym2Vec(EMtree[[m]]$probaCond),
      xlim = c(0, 1),
      ylim = c(0, 1),
      xlab = "marginales",
      ylab = "conditionnelles",
      main = m
    )
    abline(0, 1)
  })
}
load(paste0(data.dir, data.name, '-EMtree.Rdata'))

dataL <-
  data.frame(
    L = EMtree[[1]]$L,
    modele = rep(1, length(EMtree[[1]]$L)),
    index = 1:length(EMtree[[1]]$L)
  )
sapply(2:M, function(m) {
  dataL <<- rbind(dataL,
                  data.frame(
                    L = EMtree[[m]]$L,
                    modele = rep(m, length(EMtree[[m]]$L)),
                    index = 1:length(EMtree[[m]]$L)
                  ))
})
dataL %>% 
  filter(modele==4) %>% 
  ggplot(aes(index, L, color = as.factor(modele))) + theme_minimal() +
  geom_point() + geom_line() +
  coord_cartesian(xlim = c(0, 40)) +  scale_color_nejm()

# Compare edge probabilities
invisible(sapply(1:M, function(m) {
  cat(VEM[[m]]$loglik, EMtree[[m]]$L[length(EMtree[[m]]$L)], '\n')
}))
Pedge = cbind(
  F_Sym2Vec(EMtree[[1]]$probaCond),
  F_Sym2Vec(EMtree[[2]]$probaCond),
  F_Sym2Vec(EMtree[[3]]$probaCond),
  F_Sym2Vec(EMtree[[4]]$probaCond),
  F_Sym2Vec(EMtree[[5]]$probaCond)
)
colnames(Pedge) = c('M.null', 'M.date', 'M.site', 'M.all','M.inter')
par(mfrow = c(3, 2))
for (m1 in (1:(M - 1))) {
  for (m2 in ((m1 + 1):M)) {
    plot(qlogis(Pedge[, m1]), qlogis(Pedge[, m2]), xlab=m1,ylab=m2)
  }
}
Pedge %>% 
  apply(.,2,function(x) sort(x)) %>% 
  as_tibble() %>% 
  rowid_to_column() %>% 
  gather(model,proba,-rowid) %>% 
  ggplot(aes(rowid,proba,color=model))+theme_minimal()+
  geom_point() +geom_line()+ scale_color_d3()+
  coord_cartesian(xlim=c(400,528))
cor(Pedge)

# Resampling
if(Tree.res){
  Stab.seltimes = vector("list",4)
  names(Stab.seltimes)<-c("null","date","site","date_site")
  X = list(); 
  X[[1]] = matrix(1, n, 1); 
  #X[[2]] = as.matrix(lm(Data$count ~ Data$covariates$date, x=T)$x)
  X[[2]] = as.matrix(lm(Data$count ~ Data$covariates$site, x=T)$x)
  # X[[4]] = as.matrix(lm(Data$count ~ Data$covariates$date+Data$covariates$site, x=T)$x)
  #  X[[5]] = as.matrix(lm(Data$count ~ Data$covariates$date*Data$covariates$site, x=T)$x)
  T1<-Sys.time()
  #diagonale de sigma pour null et site
  invisible(sapply(1:4, function(m){
    Stab.sel3[[m]] <<- F_ResampleTreePLN(Data$count, X[[m]], matrix(0, n, p), B=500, maxIter=300,
                                         cond.tol=1e-8, cores=1)
  }))
  T2<-Sys.time()
  difftime(T2,T1)
  save(Stab.sel3, file = paste0(data.dir, data.name, '-Stab.sel_times.Rdata'))
}
unlist(lapply(Stab.sel3,median)) # temps de simu all models to full convergence
do.call(rbind,lapply(Stab.sel5, function(x){  summary(as.numeric(x$times))}))
do.call(rbind,lapply(Stab.sel5, function(x){  summary(as.numeric(x$maxIter))}))


times=gather(data.frame(t(do.call(rbind,lapply(Stab.sel5, function(x){  x$times})))),model,times)
maxIter=gather(data.frame(t(do.call(rbind,lapply(Stab.sel5, function(x){  x$maxIter})))),model,maxIter)
data<-as_tibble(cbind(times,maxIter=maxIter$maxIter))
data$normalized=data$times/data$maxIter
p1<- ggplot(data,aes(x=model,y=as.numeric(times),color=model))+
  geom_boxplot()
p2<- ggplot(data,aes(x=model,y=maxIter,color=model))+
  geom_boxplot()
p3<-ggplot(data,aes(x=model,y=normalized,color=model))+
  geom_boxplot()

grid.arrange(p1,p2,p3,nrow=3)
################################################################################
################################################################################
################################################################################

load(paste0(data.dir, data.name, '-StabSel.Rdata'))
# dans stab.sel : Pmat, alpha et iter.
Pmat_cv<-lapply(Stab.sel, function(x) {x[["Pmat"]]})
# Edge selection and comparisons
par(mfrow = c(2, 2))
edge.sel = freq.sel = list()
freq.sel<-matrix(0,nrow=ncol(Pmat_cv), ncol=M)
invisible(sapply(1:M, function(m) {
  freq.sel[,m] <<- colMeans(1 * (Pmat_cv[[m]] > 2 / p))
  edge.sel[,m] <<- 1 * (freq.sel[[m]] > freq.sel.thres)
}))
invisible(sapply(1:M, function(m) {
  if (m == 1) {
    plot(sort(freq.sel[[m]]), type = 'b')
    abline(h = freq.sel.thres)
  } else{
    points(sort(freq.sel[[m]]), col = m, type = 'b')
  }
}))

par(mfrow = c(2, 2))
invisible(sapply(1:(M - 1), function(m1) {
  sapply((m1 + 1):M, function(m2) {
    print(table(edge.sel[[m1]], edge.sel[[m2]]))
    edge.col = 1 + (2 * edge.sel[[m1]] + edge.sel[[m2]])
    plot(qlogis(Pedge[, m1]), qlogis(Pedge[, m2]), col = edge.col)
    abline(h = qlogis(2 / p), v = qlogis(2 / p))
  })
}))


# Networks
node.coord = gplot(F_Vec2Sym(edge.sel[[M]]), gmode = 'graph')
par(mfrow = c(2, 2))
invisible(sapply(1:M, function(m) {
  gplot(
    F_Vec2Sym(edge.sel[[m]]),
    gmode = 'graph',
    coord = node.coord,
    label = species,
    main = sum(edge.sel[[m]])
  )
}))

########################
# post treatment of Stability selection
##########
load(paste0(data.dir, data.name, '-StabSel.Rdata'))
load(paste0(data.dir, data.name, '-StabSel1.Rdata'))
load(paste0(data.dir, data.name, '-StabSel3.Rdata'))

Fatlist=list(Stab.sel,Stab.sel1,Stab.sel3)
names(Fatlist)<-c("cv","1","3")
# Fatlist2<-Fatlist %>% map2_dfr(.,names(.),function(x,y){
#   names(x) <- c("model1","model2","model3","model4")
#   as_tibble(x) %>% mutate(nom = y) %>% rownames_to_column() %>% 
#     gather(key, value, -rowname, -nom) %>% 
#     spread(rowname,value ) %>% mutate(`1` = map(`1`,~as.data.frame(.)))
# }) %>% unnest()%>% 
#   gather(key = "B", value = "Pmat", -key,-nom, -`2`,  -`3`)
# 
# Fatlist3<-Fatlist2 %>% 
#   rename(depth=nom,model=key, alpha=`2`, iter=`3`) %>% 
#   group_by(depth, model) 
# 
# Pmat<-lapply(Fatlist, function(x) {x[["Pmat"]]})
# 
# Fatlist2<-Fatlist %>% map2_dfr(.,names(.),function(x,y){
#   names(x) <- c("model1","model2","model3","model4")
#   as_tibble(x) %>% mutate(nom = y) %>% rownames_to_column() %>% 
#     gather(key, value, -rowname, -nom) %>% 
#     spread(rowname,value ) %>% mutate(`1` = map(`1`,function(x){
#       colMeans(1 * (x> 2 / p))
#     }))
# }) %>% unnest()%>% 
#   gather(key = "B", value = "Pmat", -key,-nom, -`2`,  -`3`)

########################
# tableau fre.sel global
#########################
build_freq<-function(list){
  df<-data.frame(freq=double(), model=character(), depth=character(), index_edge=integer())
  lapply(seq_along(list), function(j){
    depth<-names(list)[j]
    stabsel<-list[[j]]
    names(stabsel)<-c('~ 1', '~ date', '~ site', '~ date + site')
    lapply(seq_along(stabsel), function(i) { 
      #browser()
      model<-names(stabsel)[i]
      mat<-stabsel[[i]][[1]] # mat B*nb(edges)
      freq<-colMeans(1 * (mat> 2 / p)) # freq nb(edges)
      df<<-rbind(df,data.frame(freq=freq, model=model, 
                               depth=depth,index_edge=seq_along(freq)))
      
    })
  }) 
  return(df)
}
Fatfreq<-as_tibble(build_freq(Fatlist))
saveRDS(Fatfreq,"Fatfreq.rds")
Fatfreq<-readRDS("Fatfreq.rds")

DiffFreq<-Fatfreq %>% 
  spread(depth,freq) %>% 
  mutate(diff1=cv-`1`,diff3=cv-`3`, diffstart=`1`-`3`) %>% 
  gather(compare, diffs, -model,-index_edge,-cv,-`1`,-`3`) 


DiffFreq %>% 
  gather(depth, freq, -model,-index_edge, -compare, -diffs,-cv) %>%
  mutate(color=as.factor(1*(cv>=freq.sel.thres)) )%>% 
  ggplot( aes(cv, freq))+
  geom_point(aes(color=depth, pch=depth))+
  theme_minimal()+
  geom_abline()+
  geom_hline(yintercept = freq.sel.thres, size=0.5, linetype="dashed")+
  geom_vline(xintercept = freq.sel.thres, size=0.5, linetype="dashed")+
  scale_color_viridis("Iteration:", discrete=TRUE,alpha = 1, begin = 0.2, end = 0.72, direction = 1,
                      option = "A")+
  scale_shape_manual("Iteration:",values=c(19,15))+
  scale_x_continuous(breaks = c(0, 0.25,0.5,0.75,1), 
                     labels = c("0.0","0.25","0.5","0.75","1.0"))+
  facet_wrap(~model,ncol=4)+
  labs(x="Selection frequency after convergence",y="Selection frequency")+
  theme(strip.text = element_text(size=11))+ theme(panel.spacing = unit(1, "lines"))
ggsave("variation_freqSel.png", width=8.5, height=3)



freqs<-Fatfreq %>% arrange(index_edge) %>% spread(depth,freq) %>% 
  mutate_if(is.numeric,funs(1*(.>0.8))) 

freqi<-filter(freqs,model==unique(freqs$model)[4])
table(first=freqi$`3`,cv=freqi$cv) 
#####################
# Networks
############
library(ggraph)
library(tidygraph)
M<-4
p<-ncol(Data$count)
barans_net<-function(freq.sel.thres,seed){
  mat<-data.frame(F_Vec2Sym( 1*(colMeans( 1*(Stab.sel[[1]]$Pmat>2/p))>freq.sel.thres)))
  set.seed(seed)
  allNets<-tibble(P = list(mat), seuil =c("1","2","3","4") )  %>% 
    mutate(P=map( seq_along(P), function(x) {
      df<-F_Vec2Sym( 1*(colMeans( 1*(Stab.sel[[x]]$Pmat>2/p))>freq.sel.thres))
      df[lower.tri(df, diag = TRUE)]<-0
      df<-data.frame(df)
      colnames(df)<-1:ncol(df)
      df
    })) %>% 
    mutate(P = map(P,~rownames_to_column(.) %>% 
                     gather(key, value , -rowname))) %>% 
    unnest() 
  allNets<-allNets[,c(3,2,1,4)]
  groupes<-c(0,1,1,0,2,0,0,2,0,0,0,0,0,0,0,0,0,3,3,0,0,4,4,0,0,4,0,0,0,4,0,0,0)
  set_graph_style()
  
  mods<-c("Null","Date","Site","Date + Site")
  spliT<-data.frame(allNets) %>% 
    split(allNets$seuil) %>% 
    tibble(P=map(.,function(x){
      model<-switch(x$seuil[1],"1"="Null","2"="Date","3"="Site","4"="Date + Site")
      
      res<- as_tbl_graph(x, directed=FALSE) %>%
        activate(edges) %>%
        filter(value !=0) %>% 
        activate(nodes) %>% 
        mutate(group=as.factor(groupes), importance=centrality_degree(),
               keyplayer = node_is_keyplayer(k=3), model=model) 
      res %>% activate(edges) %>% 
        mutate(neibs=edge_is_incident(which(.N()$keyplayer)), model=model) %>% 
        activate(nodes) %>% 
        mutate(label=ifelse(keyplayer,name,"")) #%>% 
      #  filter(importance!=0)
    }))
  pal <- viridisLite::viridis(5, option = "C")
  pal2<-c("gray15","goldenrod1")
  lay<-create_layout(spliT$P[2][[1]],layout="circle")
  
  
  spliT$P[1][[1]] %>%
    bind_graphs(spliT$P[2][[1]] )%>%
    bind_graphs(spliT$P[3][[1]] )%>%
    bind_graphs(spliT$P[4][[1]] )%>%
    activate(nodes) %>% 
    mutate(model=factor(model,levels=mods),x=rep(lay$x,4),y=rep(lay$y,4)) %>% 
    ggraph(layout="auto")+
    geom_edge_arc(aes(color=model),curvature=0.3,show.legend=FALSE)+ 
    # scale_edge_alpha_manual(values=c(0.5,1))+
    geom_node_point(aes(color=keyplayer, size=keyplayer), show.legend=FALSE)+
    scale_edge_colour_manual(values=pal[c(1,3,2,4)], labels=mods)+
    scale_color_manual(values=pal2)+
    scale_size_manual(values=c(1.5,6))+
    geom_node_text(aes(label = label),color="black")+
    facet_nodes(~model, scales="free",ncol=4)+
    th_foreground(border=FALSE)+
    theme(strip.background = element_rect(fill="white",color="white"),
          strip.text = element_text(color="black",size=12))
}
freq.sel.thres<-0.8
seed<-200
barans_net(0.8,200)

ggsave("facet_models.png")

spliT$P[4][[1]] %>% 
  ggraph(layout="linear")+
  geom_edge_arc(aes(color=neibs),show.legend=FALSE)+ 
  geom_node_point(aes(color=keyplayer, size=keyplayer), show.legend=FALSE)+
  scale_edge_colour_manual(values=pal)+
  scale_color_manual(values=pal)+
  geom_node_text(aes(label = label),color="white")
####################
pal<-c("#196874","#a38e75","#b43434","#f0dede")
library(RColorBrewer)
spliT<-data.frame(allNets) %>% 
  split(allNets$seuil) %>% 
  tibble(P=map(.,function(x){
    model<-switch(x$seuil[1],"1"="Null","2"="Date","3"="Site","4"="Date + Site")
    
    res<- as_tbl_graph(x, directed=FALSE) %>%
      activate(edges) %>%
      filter(value !=0) %>% 
      activate(nodes) %>% 
      mutate(group=as.factor(groupes), importance=centrality_degree(),
             keyplayer = node_is_keyplayer(k=3), model=model) 
    res %>% activate(edges) %>% 
      mutate(neibs=edge_is_incident(which(.N()$keyplayer)), model=x$seuil[1]) %>% 
      activate(nodes) %>% 
      mutate(label=ifelse(keyplayer,name,"")) 
  }))
pal <- viridisLite::viridis(5, option = "C")
cols<-brewer.pal(length(mods),"OrRd")
ordCols<-cols[order(cols,(mods))]
lay<-create_layout(spliT$P[4][[1]],"fr")
spliT$P[1][[1]] %>%
  graph_join(spliT$P[2][[1]],by=c("name","group"))%>%
  graph_join(spliT$P[3][[1]],by=c("name","group"))%>%
  graph_join(spliT$P[4][[1]],by=c("name","group") ) %>% 
  activate(nodes) %>% 
  mutate(x=lay$x,y=lay$y) %>% 
  ggraph(layout="auto")+
  geom_edge_fan(aes(color=model, edge_width=model , alpha=model), spread=1.5,check_overlap = TRUE)+ 
  scale_edge_width_manual("",values=c(0.5,0.5,0.5,0.8), labels=mods)+
  scale_edge_alpha_manual("",values=c(0.5,0.4,0.8,1))+
  geom_node_point()+
  #scale_edge_colour_manual(values=brewer.pal(4,"Dark2")[c(1,3,4,2)], labels=mods)+
  scale_edge_colour_manual("",values=pal[c(4,2,3,1)], labels=mods)+
  scale_color_manual(values=pal)+
  geom_node_text(aes(label = label.y.y), nudge_y = 0.5)+
  guides( edge_alpha=FALSE)
ggsave("hairballGraph.png")

#####
#only 2 models
lay<-create_layout(spliT$P[4][[1]],"kk")
spliT$P[1][[1]] %>%
  graph_join(spliT$P[3][[1]],by=c("name","group"))%>%
  activate(nodes) %>% 
  mutate(x=lay$x,y=lay$y) %>% 
  ggraph(layout="auto")+
  geom_edge_fan(aes(color=model,alpha=..index..), check_overlap = TRUE)+ 
  geom_node_point()+
  scale_edge_colour_manual(values=brewer.pal(2,"Dark2"), labels=mods)+
  #scale_edge_colour_manual("model",values=pal[c(4,1,2,3)], labels=mods)+
  scale_color_manual(values=pal)+
  geom_node_text(aes(label = label.y), nudge_y = 0.5)+
  guides( edge_alpha=FALSE)


################
# compare tables accross models

null<-1*(colMeans( 1*(Stab.sel[[1]]$Pmat>2/p))>freq.sel.thres)
date<-1*(colMeans( 1*(Stab.sel[[2]]$Pmat>2/p))>freq.sel.thres)
site<-1*(colMeans( 1*(Stab.sel[[3]]$Pmat>2/p))>freq.sel.thres)
date_Site<-1*(colMeans( 1*(Stab.sel[[4]]$Pmat>2/p))>freq.sel.thres)

table(null,date)
table(null,site)
table(null,date_Site)
table(date,site,null) # 8 en commun aux trois
table(site, date_Site)



summary(Stab.sel[[1]]$iter)
summary(Stab.sel[[2]]$iter)
summary(Stab.sel[[3]]$iter)
summary(Stab.sel[[4]]$iter)
