% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GGMtest.R
\name{GGMtest}
\alias{GGMtest}
\title{GGMtest}
\usage{
GGMtest(
  data = X,
  edges = S,
  null_hyp = 0,
  alpha = 0.05,
  nbootstrap = 500,
  nuisance_estimaton = "lasso",
  method = "partialling out",
  DML_method = "DML2",
  s = 1,
  exponent = 1,
  penalty = list(c = 1.1),
  k_fold = 1,
  root_range = c(-100, 100)
)
}
\arguments{
\item{data}{Dataset: either matrix or dataframe}

\item{edges}{Matrix of edges for testing: each row specifies an edge}

\item{null_hyp}{A vector of null hypothesis values for the coefficients: Default is 0 for conditional independence}

\item{alpha}{Provides a value for the level of the test}

\item{nbootstrap}{Number of repetitions for the bootstrap procedure}

\item{nuisance_estimaton}{Method for nuisance parameter estimation from 'lasso', 'post-lasso' or 'sqrt-lasso'}

\item{method}{Method for point estimation, either 'root' or 'partialling out'}

\item{DML_method}{Method for point estimation, either 'DML2' or 'DML1'}

\item{s}{Number of variables combined for the confidence interval. Default is s = 1.}

\item{exponent}{Exponent for the confidence interval. Default is exponent = 1.}

\item{penalty}{Additional coefficient for the penalty term. Default value is c = 1.1.}

\item{k_fold}{Parameter for K-fold estimation. Default is k_fold = 1.}

\item{root_range}{Parameter for range of the root search (only relevant for method = 'root'). Default is root_range = (-100,100).}
}
\value{
A list with components
\item{estimates}{A vector of point estimates.}
\item{edge_list}{The matrix containing the corresponding edges (equal to input).}
\item{pvalue_max}{P-value of the hypothesis.}
\item{hyp_max}{Maximum statistic of the hypothesis.}
}
\description{
Testing conditional independence hypothesis for a gaussian graphical model.
}
\examples{
library("huge")
library("igraph")
library("GGMtest")

set.seed(42)

# generate data (different graph structures: "random", "hub", "cluster", "band" and "scale-free")
L <- huge.generator(n = 100, d = 10, graph = "cluster", g = 4)

# true Graph
true_graph <- graph_from_adjacency_matrix(as.matrix(L$theta), mode='undirected', diag=FALSE)
plot(true_graph, usearrows = FALSE, label=1:10, displaylabels=TRUE, main = "True Graph",layout= layout.fruchterman.reingold, edge.width = 2, edge.color = "black")

# index pairs for inference
S <- matrix(c(1,2,2,3,4,5), byrow = TRUE, ncol = 2)

# perform test
ggm_model <- GGMtest(data = L$data,edges = S,alpha = 0.05,nuisance_estimaton = "lasso")

# p-value:
ggm_model$pvalue_max

# plot the confidence intervals (on a subset of edges)
plot_GGMtest(ggm_model,edges = S)

}
\seealso{
\code{\link{confint.GGMtest}} for confidence intervals, \code{\link{plot_GGMtest}} for plotting options
 and \code{\link{adj_GGMtest}} for the adjacency matrix
}
