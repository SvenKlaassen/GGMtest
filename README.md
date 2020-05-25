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

As a small example simulate an example with `n = 100` observations of `p = 10` dimensional multivariate normal distribution.

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
edgeset_1 <- matrix(c(1,2,2,3,2,4,1,4,4,5), byrow = T, ncol = 2)

# perform test:
ggm_model <- GGMtest(data = L$data,
                     edges = edgeset_1,
                     alpha = 0.05,
                     nuisance_estimaton = "lasso")

# p-value:                     
ggm_model$pvalue_max

# plot the confidence intervals (on a subset of edges)
plot_GGMtest(ggm_model,edges = edgeset_1)
```

![alt text](https://github.com/SvenKlaassen/GGMtest/blob/master/plots/confidence_intervals_1.png "Confidence Intervals")

Next lets take a look at a more complex example by increasing `n = 2500` and `p = 25`. This time we will add all egdes for inference and plot the network.
```R
L <- huge.generator(n = n, d = p, graph = "random")

# all indices for inference
edgeset_2 <- t(combn(p,2))
ggm_model_2 <- GGMtest(L$data,edgeset_2, alpha = 0.05,nuisance_estimaton ="lasso")

# true Graph
par(mfrow=c(1,3))
graph_layout <- layout.circle 
#graph_layout <- layout.fruchterman.reingold
true_adj_matrix <- as.matrix(L$theta)
true_graph <- graph_from_adjacency_matrix(true_adj_matrix , mode='undirected', diag=F )
plot(true_graph, usearrows = FALSE, vertex.label=1:p,displaylabels=T,main = "True Graph",
     layout= graph_layout,edge.color = "black",edge.width = 2)

# estimated graph
est_adj_matrix <- adj_GGMtest(ggm_model_2)
est_graph <- igraph::graph_from_adjacency_matrix(est_adj_matrix , mode='undirected', diag=F )
plot(est_graph, usearrows = FALSE, vertex.label=1:p,displaylabels=T,main = "Estimated Graph",
     layout= graph_layout,edge.color = "black",edge.width = 2)

#calculation matrix difference
diff_matrix <- true_adj_matrix - est_adj_matrix 
diff_graph <- igraph::graph_from_adjacency_matrix(diff_matrix  , mode='undirected', diag=F, weighted = T)
E(diff_graph)$color[E(diff_graph)$weight == 1] <- "blue"
E(diff_graph)$color[E(diff_graph)$weight == -1] <- "red"
plot(diff_graph, usearrows = FALSE, vertex.label=1:p,displaylabels=T,main = "Graph Differences",
     layout= graph_layout,edge.width = 2)
legend(x=-1.5, y=-1.1, c("False Positive","False Negative"), pch=21,
       col="#777777", pt.bg=c("red","blue"), pt.cex=2, cex=.8, bty="n", ncol=1)

par(mfrow=c(1,1))
```
![alt text](https://github.com/SvenKlaassen/GGMtest/blob/master/plots/example_large_1.png "Estimated GGM")

## Application
A small application from psychology. The  `bfi` dataset contains `n = 2800` oberservations of 25 personality self report intems. For a documentation on the covariables see `bfi` in the `psych` package.
```R
library("psych")
library("dplyr")

df <- bfi %>%
  select(1:25) %>%
  na.omit() %>%
  apply(2,function(x) x-mean(x))

# all indices for inference
edgeset_3 <- t(combn(25,2))

ggm_model_3 <- GGMtest(df,edgeset_3,alpha = 0.05,nuisance_estimaton ="lasso")

# estimated graph
est_adj_matrix <- adj_GGMtest(ggm_model_3)
est_graph <- igraph::graph_from_adjacency_matrix(est_adj_matrix , mode='undirected', diag=F )
plot(est_graph, usearrows = FALSE, vertex.label=colnames(df),displaylabels=T,main = "Estimated Graph",
     layout= layout.fruchterman.reingold,edge.color = "black",edge.width = 2)

```
![alt text](https://github.com/SvenKlaassen/GGMtest/blob/master/plots/application_bfi_1.png "Estimated GGM")

## References

* Klaassen, S., Kück, J., Spindler, M., & Chernozhukov, V. (2018). Uniform inference in high-dimensional Gaussian graphical models. arXiv preprint arXiv:1808.10532.
* Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters.
* Belloni, Alexandre, et al. "Uniformly valid post-regularization confidence regions for many functional parameters in z-estimation framework." Annals of statistics 46.6B (2018): 3643.

