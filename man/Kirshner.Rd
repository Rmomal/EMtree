% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunctionsTree.R
\name{Kirshner}
\alias{Kirshner}
\title{Computing edge probabilities using Kirshner (07) formula}
\usage{
Kirshner(W)
}
\arguments{
\item{W}{Squared weight matrix}
}
\value{
Edges probabilities as defined in Kirshner 2007
}
\description{
Computing edge probabilities using Kirshner (07) formula
}
\examples{
W = matrix(c(1,1,3,1,1,1,3,1,1),3,3,byrow=TRUE)
Kirshner(W)
}
