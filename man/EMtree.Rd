% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_inference.R
\name{EMtree}
\alias{EMtree}
\title{Core computing function}
\usage{
EMtree(
  PLN.Cor,
  n = NULL,
  maxIter = 30,
  unlinked = NULL,
  random.init = FALSE,
  cond.tol = 1e-10,
  eps = 0.001,
  verbatim = TRUE,
  plot = FALSE
)
}
\arguments{
\item{PLN.Cor}{Either an object resulting from the use of the `PLN` function from package `PLNmodels`, or an estimation of Gaussian data correlation matrix.}

\item{n}{Number of samples, required if a correlation matrix is supplied}

\item{maxIter}{Maximum number of iterations for EMtree}

\item{unlinked}{An optional vector of nodes which are not linked with each other}

\item{random.init}{A boolean for trying a random initialization of the EM}

\item{cond.tol}{Tolerance parameter for the conditioning of psi matrix}

\item{eps}{Precision parameter controlling the convergence of weights beta}

\item{verbatim}{Talks if set to TRUE}

\item{plot}{Plots likelihood if set to TRUE}
}
\value{
\itemize{
 \item{edges_prob: }{p x p matrix of edges probabilities}
 \item{edges_weight: }{p x p matrix of edges weights for any spanning tree}
 \item{logpY: }{vector of log-likelihoods}
 \item{maxIter: }{final number of iterations EMtree has ran}
 \item{timeEM: }{EMtree computation time}
 \item{alpha: }{data signal/noise ratio}
}
}
\description{
Core computing function
}
\examples{
n=30
p=10
Y=data_from_scratch("tree",p=p)$data
PLN_Y = PLNmodels::PLN(Y~1)
EMtree(PLN.Cor=PLN_Y,verbatim=TRUE, plot=TRUE)
}
