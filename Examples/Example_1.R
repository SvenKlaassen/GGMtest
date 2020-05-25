rm(list=ls())

package_list <- c("huge","tictoc","igraph","GGMtest")
invisible(lapply(package_list, library, character.only = TRUE))


set.seed(42)
#### generate data ####
n <- 300  # number of observations
p <- 28   # dimension of the gaussian vector
l <- 200    # number of independent Monte-Carlo Estimation


# generate data
L <- huge.generator(n = l*n, d = p, graph = "cluster", g = 4) #different graph structures including "random", "hub", "cluster", "band" and "scale-free".
# true Graph
true_graph <- graph_from_adjacency_matrix(as.matrix(L$theta) , mode='undirected', diag=F )
plot(true_graph, usearrows = FALSE, label=1:p,displaylabels=T,main = "True Graph",
     layout= layout.fruchterman.reingold,edge.width = 2,edge.color = "black")

# index pairs for inference
# Testing independence of cluster one and two
p_2 = p/4
S_1 <- c()
for (i in 1:p_2){
  S_1 <- c(S_1,rep(i,p_2))
}
for (i in 1:p_2){
  S_1 <- c(S_1,(p_2+1):(2*p_2))
}
S <- matrix(S_1,byrow =F,ncol =2)
p_1 = dim(S)[1]

#### Estimation ####
v_hyp_max <- vector("logical",l)

tic("total time")
tic("time difference")
for (h in 1:l){
  X <- L$data[((h-1)*n+1):(h*n),]

  #### GGMtest ####
  ggm_model <- GGMtest(data = X,
                       edges = S,
                       null_hyp = rep(0,dim(S)[1]),
                       alpha = 0.1,
                       nbootstrap = 500,
                       nuisance_estimaton = "lasso",
                       s=1,
                       exponent=1,
                       k_fold = 1,
                       rnd_seed = NULL)

  v_hyp_max[h]<- ggm_model$hyp_max

  if (h%%(l/10)==0){
    cat("step", h, "out of" , l,"\n")
    toc()
    tic("time difference")
  }
}
toc(quiet = T )
toc()

sprintf("Coverage: %s", sum(v_hyp_max)/l)
