% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunctionsTree.R
\name{Laplacian}
\alias{Laplacian}
\title{Compute a matrix Laplacian}
\usage{
Laplacian(W)
}
\arguments{
\item{W}{Squared weight matrix}
}
\value{
The matrix Laplacian of W
}
\description{
Compute a matrix Laplacian
}
\examples{
uniformeW=matrix(1,3,3)
Laplacian(uniformeW)
}
