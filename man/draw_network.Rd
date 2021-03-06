% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_Graphs.R
\name{draw_network}
\alias{draw_network}
\title{Plots a network}
\usage{
draw_network(
  adj_matrix,
  title = "",
  label_size = 4,
  curv = 0,
  width = 1,
  shade = FALSE,
  remove_isolated = FALSE,
  btw_rank = 2,
  layout = NULL,
  stored_layout = NULL,
  nodes_label = NULL,
  nodes_size = c(2, 5),
  pal_edges = NULL,
  pal_nodes = NULL,
  node_groups = NULL,
  edge_groups = NULL,
  legend = FALSE,
  leg.pos = "bottom",
  title.family = "sans",
  title.size = 12
)
}
\arguments{
\item{adj_matrix}{graph adjacency matrix}

\item{title}{graph title}

\item{label_size}{size of labels}

\item{curv}{edges curvature}

\item{width}{maximum width for the edges}

\item{shade}{if TRUE, shades the edges unlinked to nodes with high betweenness}

\item{remove_isolated}{Removes isolated nodes if TRUE}

\item{btw_rank}{betweenness rank +1 of highlighted nodes. If set to 1, none are highlighted.}

\item{layout}{optional layout from `ggraph` or `ìgraph` package.}

\item{stored_layout}{optional data.frame of point positions from a previous graph, possibly created using ggraph::create_layout()}

\item{nodes_label}{optional labels for nodes.}

\item{nodes_size}{size of nodes, possibility to specify a size per group}

\item{pal_edges}{optional palette for edges}

\item{pal_nodes}{optional palette for nodes}

\item{node_groups}{optional vector separating the nodes into groups}

\item{edge_groups}{optional matrix specifying groups of edges}

\item{legend}{optional boolean for generating a legend when groups are provided}

\item{leg.pos}{optional legend position ("bottom", "top", "right" or "left")}

\item{title.family}{font family of the title}

\item{title.size}{font size of the title}
}
\value{
\itemize{
\item{G} {the network as a ggplot2 object, with highlighted high betweenness nodes}
\item{graph_data}{ data needed for plotting the network}
\item{legend}{ ggplot2 information about the legend}
}
}
\description{
Plots a network
}
\examples{
set.seed(1)
adj_matrix= SimCluster(p=20,k=2,dens=0.5, r=10)
groups=sample(3,20, replace=TRUE)
draw_network(adj_matrix,"Cluster graph", layout="stress",
shade=TRUE, btw_rank=3)$G
draw_network(adj_matrix,"Cluster with groups of nodes",
layout="stress", shade=TRUE, btw_rank=3,
node_groups=groups, legend=TRUE,title.family="mono",title.size=12)$G
}
