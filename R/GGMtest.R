#' GGMtest
#'
#' Testing conditional independence hypothesis for a gaussian graphical model.
#'
#' @param data Dataset: either matrix or dataframe
#' @param edges Matrix of edges for testing: each row specifies an edge
#' @param null_hyp A vector of null hypothesis values for the coefficients: Default is 0 for conditional independence
#' @param alpha Provides a value for the level of the test
#' @param nbootstrap Number of repetitions for the bootstrap procedure
#' @param nuisance_estimaton Method for nuisance parameter estimation from 'lasso', 'post-lasso' or 'sqrt-lasso'
#' @param method Method for point estimation, either 'root' or 'partialling out'
#' @param DML_method Method for point estimation, either 'DML2' or 'DML1'
#' @param s Number of variables combined for the confidence interval. Default is s = 1.
#' @param exponent Exponent for the confidence interval. Default is exponent = 1.
#' @param penalty Additional coefficient for the penalty term. Default value is c = 1.1.
#' @param k_fold Parameter for K-fold estimation. Default is k_fold = 1.
#' @param root_range Parameter for range of the root search (only relevant for method = 'root'). Default is root_range = (-100,100).
#'
#' @return A list with components
#' \item{estimates}{A vector of point estimates.}
#' \item{edge_list}{The matrix containing the corresponding edges (equal to input).}
#' \item{pvalue_max}{P-value of the hypothesis.}
#' \item{hyp_max}{Maximum statistic of the hypothesis.}
#'
#' @examples
#' library("huge")
#' library("igraph")
#' library("GGMtest")
#'
#' set.seed(42)
#'
#' # generate data (different graph structures: "random", "hub", "cluster", "band" and "scale-free")
#' L <- huge.generator(n = 100, d = 10, graph = "cluster", g = 4)
#'
#' # true Graph
#' true_graph <- graph_from_adjacency_matrix(as.matrix(L$theta), mode='undirected', diag=FALSE)
#' plot(true_graph, usearrows = FALSE, label=1:10, displaylabels=TRUE, main = "True Graph",layout= layout.fruchterman.reingold, edge.width = 2, edge.color = "black")
#'
#' # index pairs for inference
#' S <- matrix(c(1,2,2,3,4,5), byrow = TRUE, ncol = 2)
#'
#' # perform test
#' ggm_model <- GGMtest(data = L$data,edges = S,alpha = 0.05,nuisance_estimaton = "lasso")
#'
#' # p-value:
#' ggm_model$pvalue_max
#'
#' # plot the confidence intervals (on a subset of edges)
#' plot_GGMtest(ggm_model,edges = S)
#'
#' @seealso \code{\link{confint.GGMtest}} for confidence intervals, \code{\link{plot_GGMtest}} for plotting options
#'  and \code{\link{adj_GGMtest}} for the adjacency matrix
#'
#' @export
#'
#'

GGMtest <- function(data = X,
                    edges = S,
                    null_hyp = 0,
                    alpha = 0.05,
                    nbootstrap = 500,
                    nuisance_estimaton = 'lasso',
                    method = 'partialling out',
                    DML_method = 'DML2',
                    s = 1,
                    exponent = 1,
                    penalty = list(c = 1.1),
                    k_fold = 1,
                    root_range = c(-100,100)) {
  X <- as.matrix(data)
  S <- as.matrix(edges)
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (length(null_hyp) == 1){
    beta_0 <- rep(null_hyp,dim(S)[1])
  } else {
    beta_0 <- null_hyp
  }


  #### Checking Arguments ####

  checkmate::checkNumeric(alpha,0,1)
  checkmate::checkChoice(nuisance_estimaton, c("lasso","post-lasso","sqrt-lasso"))
  checkmate::checkChoice(method, c('partialling out','root'))

  if (dim(S)[2] != 2){
    stop("Invalid argument: S has to be a matrix with 2 columns.")
  }
  if (any(S[,1]==S[,2])){
    stop("Invalid argument: Not a valid edge. The indices cant be identical.")
  }
  if (!isTRUE(all.equal(as.vector(S), as.integer(S)))){
    stop("Invalid argument: The indices have to be integers.")
  }
  if (max(as.integer(S)>p)){
    stop("Invalid argument: At least one index is larger than the number of regressors.")
  }
  if (DML_method == 'DML2' & method == 'root'){
    stop("Only partialling out is implemented for DML2.")
  }


  #### K-fold partition ####
  num_drop <- n-k_fold*floor(n/k_fold) # number of indices to be added seperately

  if (num_drop > 0) {
    indices_drop <- sample(n,num_drop) # random indices to be added seperately
    indices_remain <- (1:n)[-indices_drop] # the remaining indices
    folds <- split(sample(indices_remain,n-num_drop,replace=FALSE), as.factor(1:k_fold)) # calculation of the preversion of the folds

    # adding the dropped indices
    for (i in 1:length(indices_drop)){
      folds[[i]] <- c(folds[[i]],indices_drop[i])
    }
  } else {
    folds <- split(sample(n,n,replace=FALSE), as.factor(1:k_fold)) # calculation of the preversion of the folds
  }
  #### General variables and score function: ####

  p1 <- dim(S)[1]  # number of estimations


  block1 <- unique(S[,1])

  Score <- function(theta,X,r1,r2){
    score <- (r1-theta*X)*r2
    return(score)
  }

  beta_mat <- array(NA,dim=c(p1,k_fold))

  pred1 <- array(NA,dim=c(p1,(n/k_fold+1),k_fold)) # empty entries possible
  pred2 <- array(NA,dim=c(p1,(n/k_fold+1),k_fold)) # empty entries possible

  # definition subsamples

  for (num_fold in 1:k_fold){
    if (k_fold == 1){
      sample1_index <- 1:n
      sample2_index <- 1:n
    } else {
      sample2_index <- folds[[num_fold]] # indices for estimation of the score function
      sample1_index <- setdiff(1:n,folds[[num_fold]]) # indices for estimation of the nuisance function
    }
    n1 <- length(sample1_index)
    n2 <- length(sample2_index)

    gamma <- 0.1/log(n1) # because of the p1-adjustment

    lambda_1 <- penalty$c*1/(2*sqrt(n1))*stats::qnorm(1-gamma/(2*p1*max((p-1),n1))) #penalty term eta one
    lambda_2 <- penalty$c*1/(2*sqrt(n1))*stats::qnorm(1-gamma/(2*p1*max((p-2),n1))) #penalty term eta two

    #### Estimation (includes elements for variance estimation) ####

    for (i in block1){ # estimation relevant coefficients
      Y1 <- X[sample1_index,i]
      X1 <- X[sample1_index,-i]
      if (nuisance_estimaton == "lasso"){
        r1 <- hdm::rlasso(X1,Y1, post=FALSE,intercept=F,penalty = list(c = penalty$c,gamma = gamma/p1))$coefficients
      } else if (nuisance_estimaton == "post-lasso") {
        r1 <- hdm::rlasso(X1,Y1, post=TRUE,intercept=F,penalty = list(c = penalty$c,gamma = gamma/p1))$coefficients
      } else if (nuisance_estimaton == "sqrt-lasso") {
        r1 <- picasso::picasso(X1,Y1, family="sqrtlasso",lambda = lambda_1)$beta[,1]
      }
      #get the corresponding indices
      block2 <- S[,2][S[,1]== i]

      for (j in block2){

        index <- match(TRUE, i == S[,1] & j == S[,2])

        Y2 <- X[sample1_index,j]
        X2 <- X[sample1_index,-c(i,j)]

        if (nuisance_estimaton == "lasso"){
          r2 <- hdm::rlasso(X2,Y2, post=FALSE,intercept=F,penalty = list(c = penalty$c,gamma = gamma/p1))$coefficients
        } else if (nuisance_estimaton == "post-lasso") {
          r2 <- hdm::rlasso(X2,Y2, post=TRUE,intercept=F,penalty = list(c = penalty$c,gamma = gamma/p1))$coefficients
        } else if (nuisance_estimaton == "sqrt-lasso") {
          r2 <- picasso::picasso(X2,Y2, family="sqrtlasso",lambda = lambda_2)$beta[,1]
        }

        if (j > i) {
          c1 <- -1
        } else {
          c1 <- 0
        }

        pred1[index,1:n2,num_fold] <- X[sample2_index,i]- X[sample2_index,-c(i,j)]%*%r1[-(j+c1)]

        pred2[index,1:n2,num_fold] <- X[sample2_index,j]- X[sample2_index,-c(i,j)]%*%r2 #j=k

        #### Estimator ####
        if (DML_method == "DML1"){
          if (method == 'partialling out'){
            #partialling out (first order equivalent)
            beta1 <- stats::lm(pred1[index,1:n2,num_fold]~pred2[index,1:n2,num_fold])$coefficients[2]
          } else if (method == 'root'){
            beta1 <- stats::uniroot(function(x) mean(Score(x,X[sample2_index,j],pred1[index,1:n2,num_fold],pred2[index,1:n2,num_fold])),interval = c(root_range[1],root_range[2]))$root
          }
          beta_mat[index,num_fold] <- beta1
        }
      }
    }
  }
  if (DML_method == "DML1"){
    beta_vec <- apply(beta_mat,1,mean) #calculating the mean over all estimators for each component of beta
  } else if (DML_method == "DML2"){
    beta_vec <- array(NA,dim=c(p1))
    for (i in block1){
      #get the corresponding indices
      block2 <- S[,2][S[,1]== i]
      for (j in block2){
        index <- match(TRUE, i == S[,1] & j == S[,2])
        error_1 <- t(R.utils::wrap(pred1[index,1:n2,,drop = F],map=list(1,NA)))
        error_2 <- t(R.utils::wrap(pred2[index,1:n2,,drop = F],map=list(1,NA)))
        beta_vec[index] <- stats::lm(error_1 ~ error_2)$coefficients[2]
      }
    }
  }

  #### Estimation of J ####

  Gamma <- array(NA,dim=c(p1,k_fold))

  for (num_fold in 1:k_fold){
    if (k_fold == 1){
      sample1_index <- 1:n
      sample2_index <- 1:n
    } else {
      sample2_index <- folds[[num_fold]] # indices for estimation of the score function
      sample1_index <- setdiff(1:n,folds[[num_fold]]) # indices for estimation of the nuisance function
    }
    n1 <- length(sample1_index)
    n2 <- length(sample2_index)

    for (i in block1){
      block2 <- S[,2][S[,1]== i]
      for (j in block2){
        index <- match(TRUE, i == S[,1] & j == S[,2])
        Gamma[index,num_fold] <- mean((-X[sample2_index,j])*pred2[index,1:n2,num_fold])
      }
    }
  }
  Gamma_vec <- apply(Gamma,1,mean)


  #### Estimation of sigma ####

  sigma_est <- array(NA,dim=c(p1,k_fold))
  for (num_fold in 1:k_fold){
    if (k_fold == 1){
      sample1_index <- 1:n
      sample2_index <- 1:n
    } else {
      sample2_index <- folds[[num_fold]] # indices for estimation of the score function
      sample1_index <- setdiff(1:n,folds[[num_fold]]) # indices for estimation of the nuisance function
    }
    n1 <- length(sample1_index)
    n2 <- length(sample2_index)
    score_square <-array(NA,dim=c(p1,n2))

    for (i in block1){
      block2 <- S[,2][S[,1]== i]

      for (j in block2){
        index <- match(TRUE, i == S[,1] & j == S[,2])

        score_square[index,] <- (Score(beta_vec[index],X[sample2_index,j],pred1[index,1:n2,num_fold],pred2[index,1:n2,num_fold]))^2
        sigma_est[index,num_fold] <- sqrt(mean(Gamma_vec[index]^(-2)*score_square[index,]))
      }
    }
  }
  sigma_vec <- apply(sigma_est,1,mean)


  #### Estimation of psi ####

  psi_est <- array(NA,dim=c(p1,(n/k_fold+1),k_fold))


  for (num_fold in 1:k_fold){
    if (k_fold == 1){
      sample1_index <- 1:n
      sample2_index <- 1:n
    } else {
      sample2_index <- folds[[num_fold]] # indices for estimation of the score function
      sample1_index <- setdiff(1:n,folds[[num_fold]]) # indices for estimation of the nuisance function
    }
    n1 <- length(sample1_index)
    n2 <- length(sample2_index)
    score_square <-array(NA,dim=c(p1,n2)) # drop!

    for (i in block1){
      block2 <- S[,2][S[,1]== i]

      for (j in block2){
        index <- match(TRUE, i == S[,1] & j == S[,2])
        psi_est[index,1:n2,num_fold] <- -(sigma_vec[index]*Gamma_vec[index])^(-1)*Score(beta_vec[index],X[sample2_index,j],pred1[index,1:n2,num_fold],pred2[index,1:n2,num_fold])

      }
    }
  }


  #### Multiplier Bootstrap ####

  dist_est <- array(NA,dim= c(p1,nbootstrap))

  for (b in 1:nbootstrap){
    if (num_drop == 0){
      epsilon <- matrix(stats::rnorm(n+k_fold),byrow=F,ncol=k_fold)
    } else {
      epsilon <- matrix(stats::rnorm(ceiling(n/k_fold)*k_fold),byrow=F,ncol=k_fold)
    }
    dist_est[,b] <- apply(psi_est,1,function(x) n^(-1/2) * sum(epsilon * x, na.rm = T))
    #for (j in 1:p1){dist_est[j,b] <- n^(-1/2)*sum(epsilon*psi_est[j,,],na.rm = T)}
  }

  #### Hypotheses Testing ####
  stat1 <- sqrt(n)*(beta_vec-beta_0)*sigma_vec^(-1)


  #### squared over some combinations ####

  dist_est2 <- array(0, dim=c(p1,nbootstrap))
  stat2 <- array(0, dim= c(p1))
  vec1 <- sample(p1,p1)
  count_dist_est2 <- 1

  # defining the function of the sum to the power of exponent
  sum_exponent <- function(x,s,start_ind,exponent,ind_vec){
    ind_vec <- c(ind_vec,ind_vec[1:(s-1)]) # the adjust for the end of the vector
    result <- x[ind_vec[start_ind]]^exponent
    if (s > 1){
      for (u in 1:(s-1)){
        result <- result + x[ind_vec[start_ind+u]]^exponent
      }
    }
    return(result)
  }

  for (j in 1:p1){
    dist_est2[count_dist_est2,] <- apply(dist_est,2,function(x) sum_exponent(x,s,start=j,exponent=exponent,ind_vec= vec1))
    stat2[count_dist_est2]<- sum_exponent(stat1,s,start=j,exponent=exponent,ind_vec= vec1)
    count_dist_est2 = count_dist_est2+1
  }
  nsample=apply(abs(dist_est2),2,max)
  quant_est <- stats::quantile(nsample,probs = c(1-alpha,1-alpha/2,alpha/2))

  hyp_max <- (max(abs(stat2))<= quant_est[[1]])

  hyp_sphere <- (max(abs(stat2))<= quant_est[[2]] && quant_est[[3]] <= max(abs(stat2)))

  Nsample <- sort(nsample)
  target_est <- max(abs(stat2))

  #### p-value ####

  hmax=TRUE
  b=0
  while (hmax==TRUE && b<=nbootstrap){
    hmax=(target_est<= Nsample[(nbootstrap-b)])
    b=b+1

  }
  pvalue_max=(b-1)/nbootstrap

  hsphere=TRUE
  b=0
  while (hsphere==TRUE && b<=(nbootstrap/2)){
    hsphere=(target_est<= Nsample[(nbootstrap-b)] && Nsample[(b+1)] <= target_est)
    b=b+1

  }
  pvalue_sphere=2*(b-1)/nbootstrap

  #save additional parameters
  add_par <- list(sigma_est = sigma_vec, n = n, alpha = alpha,beta_0 = beta_0,s = s,exponent = exponent,
                  random_order = vec1, p = p)

  result <- list(estimates = beta_vec,edge_list = S, pvalue_max = pvalue_max, hyp_max=hyp_max, alpha = alpha,
                 nbootstrap = nbootstrap, nuisance_estimaton =  nuisance_estimaton, penalty = penalty, folds = folds,
                 Nsample = Nsample, target_est = target_est, hyp_sphere= hyp_sphere, pvalue_sphere = pvalue_sphere,
                 quantile =quant_est, additional_parameters = add_par)
  class(result) <- "GGMtest"
  return(result)
}
