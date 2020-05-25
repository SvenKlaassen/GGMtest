# GGMtest

Hypothesis testing for conditional independence in gaussian graphical models. This repository is the official implementation of [Uniform inference in high-dimensional Gaussian graphical models](https://arxiv.org/pdf/1808.10532.pdf).

## Description

Testing conditional independence in gaussian graphical models becomes a high-dimensional inference problem due to the large number of possible edges. To increase power in contrast to multiple testing adjustments, we employ the double machine learning framework.

## Installing GGMtest

To install this package in R, run the following commands:

```R
install.packages("devtools")
library(devtools)
install_github("SvenKlaassen/GGMtest")
```

## Examples

As a small example simulate an example with \R{n = 100} observations of \code{p = 10} dimensional multivariate normal distribution.

```R
rm(list=ls())

library("huge")
library("igraph")
library("GGMtest")

set.seed(42)

# generate data (different graph structures: "random", "hub", "cluster", "band" and "scale-free") 
L <- huge.generator(n = 100, d = 10, graph = "cluster", g = 4) 
# true Graph
true_graph <- graph_from_adjacency_matrix(as.matrix(L$theta), mode='undirected', diag=F)

plot(true_graph, usearrows = FALSE, label=1:10, displaylabels=T, main = "True Graph",
     layout= layout.fruchterman.reingold, edge.width = 2, edge.color = "black")

```

![alt text](https://github.com/SvenKlaassen/GGMtest/blob/master/plots/True_Graph_1.png "True Graph")

```R
# index pairs for inference
S <- matrix(c(1,2,2,3,2,4,1,4,4,5), byrow = T, ncol = 2)

# perform test:
ggm_model <- GGMtest(data = L$data,
                     edges = S,
                     alpha = 0.05,
                     nuisance_estimaton = "lasso")

# p-value:                     
ggm_model$pvalue_max

# plot the confidence intervals (on a subset of edges)
plot_GGMtest(ggm_model,edges = S)
```

![alt text](https://github.com/SvenKlaassen/GGMtest/blob/master/plots/confidence_intervals_1.png "True Graph")

Next lets take a look at a more complex example by increasing \code{n = 2500} and \code{p = 25}.

## Application

## References

* Klaassen, S., Kück, J., Spindler, M., & Chernozhukov, V. (2018). Uniform inference in high-dimensional Gaussian graphical models. arXiv preprint arXiv:1808.10532.
* Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters.
* Belloni, Alexandre, et al. "Uniformly valid post-regularization confidence regions for many functional parameters in z-estimation framework." Annals of statistics 46.6B (2018): 3643.

