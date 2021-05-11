#' GGMtest2
#'
#' Testing conditional independence hypothesis for a gaussian graphical model.
#'
#' @param data Dataset: either matrix or dataframe
#' @param edges Matrix of edges for testing: each row specifies an edge
#' @param nuisance_estimaton Method for nuisance parameter estimation from 'lasso', 'post-lasso' or 'sqrt-lasso'
#' @param method Method for point estimation, either 'root' or 'partialling out'
#' @param DML_method Method for point estimation, either 'DML2' or 'DML1'
#' @param k_fold Parameter for K-fold estimation. Default is k_fold = 1.
#' @param penalty Additional coefficient for the penalty term. Default value is c = 1.1.
#'
#' @return A list with components
#' \item{estimates}{A vector of point estimates.}
#' \item{edge_list}{The matrix containing the corresponding edges (equal to input).}
#' \item{sigma_est}{Estimates of the standard deviation.}
#' \item{psi_est}{Estimates of the score vector.}
#' \item{additional_parameters}{Additional parameters.}
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
#' ggm_model <- GGMtest2(data = L$data,edges = S,nuisance_estimaton = "lasso")
#'
#' # Create Confidence Region:
#' create_CR(ggm_model)
#'
#'
#' @seealso \code{\link{confint.GGMtest}} for confidence intervals, \code{\link{plot_GGMtest}} for plotting options
#'  and \code{\link{adj_GGMtest}} for the adjacency matrix
#'
#' @export
#'
#'

GGMtest2 <- function(data = X,
                     edges = S,
                     nuisance_estimaton = 'lasso',
                     method = 'robust',
                     DML_method = 'DML2',
                     k_fold = 1,
                     penalty = list(c = 1.1)) {
  #### Checking Arguments ####
  checkmate::assertChoice(nuisance_estimaton, c("lasso","post-lasso","sqrt-lasso"))
  checkmate::assertChoice(method, c('robust','partialling out'))
  checkmate::assertChoice(DML_method, c('DML1','DML2'))
  checkmate::assertIntegerish(k_fold, lower = 1)
  checkmate::assertList(penalty)

  X <- as.matrix(data)
  S <- as.matrix(edges)
  n <- dim(X)[1]
  p <- dim(X)[2]

  checkmate::assertMatrix(S,ncols = 2)
  checkmate::assertIntegerish(S,upper = p)
  if (any(S[,1]==S[,2])){
    stop("Invalid argument: Not a valid edge. The indices cant be identical.")
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
          if (method == 'robust'){
            beta1 <- mean(pred2[index,1:n2,num_fold]*pred1[index,1:n2,num_fold])/mean(pred2[index,1:n2,num_fold]*X[sample2_index,j])
          } else if (method == 'partialling out'){
            #partialling out (first order equivalent)
            beta1 <- stats::lm(pred1[index,1:n2,num_fold] ~ pred2[index,1:n2,num_fold])$coefficients[2]
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
        if (method == 'robust'){
          if (k_fold == 1){
            X_j <- X[,j]
          } else {
            X_j <- X[Reduce(c,folds),j]
          }
          beta_vec[index] <- mean(error_2*error_1)/mean(error_2*X_j)
        } else if (method == 'partialling out'){
          beta_vec[index] <- stats::lm(error_1 ~ error_2)$coefficients[2]
        }
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

  #save additional parameters
  add_par <- list(n = n,
                  p = p,
                  p1 = p1,
                  nuisance_estimaton =  nuisance_estimaton,
                  penalty = penalty,
                  k_fold = k_fold,
                  folds = folds,
                  num_drop = num_drop)

  result <- list(estimates = beta_vec,
                 edge_list = S,
                 sigma_est = sigma_vec,
                 psi_est = psi_est,
                 additional_parameters = add_par)

  class(result) <- "GGMtest2"
  return(result)
}


#' Create Confidence regions for a `GGMtest2` object.
#'
#' @param model The `GGMtest2` object.
#' @param null_hyp A vector of null hypothesis values for the coefficients: Default is 0 for conditional independence
#' @param alpha The corresponding level.
#' @param B The number of used Bootstrap repetitions.
#' @param s The number of s-sparse combinations.
#' @param exp The corresponding exponent of the s-sparse sets.
#'
#' @return A list with the following components.
#' \item{estimates}{A vector of point estimates.}
#' \item{edge_list}{The matrix containing the corresponding edges (equal to input).}
#' \item{pvalue_max}{P-Value of the maximum statistic.}
#' \item{pvalue_sphere}{P-Value of the approx. Sphere.}
#' \item{hyp_max}{`FALSE` if the hypthesis is rejected with the maximum statistic.}
#' \item{hyp_sphere}{`FALSE` if the hypthesis is rejected with the approx. Sphere.}
#' \item{vol_max}{Volume of the maximum statistic.}
#' \item{vol_sphere}{Volume of  the approx. Sphere.}
#' \item{quantiles}{Estimates quantiles.}
#'
#' @export

create_CR <- function(model, null_hyp = 0, alpha = 0.05, B = 500, s = 1, exp = 1){
  #check arguments
  checkmate::assertClass(model,"GGMtest2")
  checkmate::assertNumeric(null_hyp, min.len = 1)
  checkmate::assertNumber(alpha,lower = 0, upper = 1)
  checkmate::assertIntegerish(B,lower = 1)
  checkmate::assertIntegerish(s,lower = 1)
  checkmate::assertIntegerish(exp, lower = 1)

  #specify parameters
  n <- model$additional_parameters$n
  p1 <- model$additional_parameters$p1
  k_fold <- model$additional_parameters$k_fold
  num_drop <- model$additional_parameters$num_drop

  if (length(null_hyp) == 1){
    beta_0 <- rep(null_hyp,p1)
  } else {
    checkmate::assertNumeric(null_hyp, len = p1)
    beta_0 <- null_hyp
  }

  #### Multiplier Bootstrap ####

  dist_est <- array(NA,dim= c(p1,B))

  for (b in 1:B){
    if (num_drop == 0){
      epsilon <- matrix(stats::rnorm(n+k_fold),byrow=F,ncol=k_fold)
    } else {
      epsilon <- matrix(stats::rnorm(ceiling(n/k_fold)*k_fold),byrow=F,ncol=k_fold)
    }
    dist_est[,b] <- apply(model$psi_est,1,function(x) n^(-1/2) * sum(epsilon * x, na.rm = T))
    #for (j in 1:p1){dist_est[j,b] <- n^(-1/2)*sum(epsilon*psi_est[j,,],na.rm = T)}
  }

  #### Hypotheses Testing ####
  stat1 <- sqrt(n)*(model$estimates-beta_0)*model$sigma_est^(-1)

  #### squared over some combinations ####
  dist_est2 <- array(0, dim=c(ceiling(p1/s),B))
  stat2 <- array(0, dim= ceiling(p1/s))
  if (s == 1){
    vec1 <- subsets <- 1:p1
  } else {
    vec1 <- sample(p1,p1) #randomly reorder equations
    #create disjunct subsets of size s
    subsets <- split(vec1,factor(rep(1:ceiling(p1/s),s))[1:p1])
  }

  count_dist_est2 <- 1
  for (subset in subsets){
    #changed to lp-Ball (added absolute value)
    dist_est2[count_dist_est2,] <- apply(dist_est,2,function(x) sum(abs(x[subset])^exp))
    stat2[count_dist_est2]<- sum(abs(stat1[subset])^exp)
    count_dist_est2 = count_dist_est2+1
  }

  nsample <- apply(abs(dist_est2),2,max)
  quant_est <- stats::quantile(nsample,probs = c(1-alpha,1-alpha/2,alpha/2))

  hyp_max <- (max(abs(stat2)) <= quant_est[[1]])
  hyp_sphere <- (max(abs(stat2)) <= quant_est[[2]] && quant_est[[3]] <= max(abs(stat2)))

  Nsample <- sort(nsample)
  target_est <- max(abs(stat2))

  #### p-value ####
  hmax=TRUE
  b=0
  while (hmax==TRUE && b<=B){
    hmax=(target_est <= Nsample[(B-b)])
    b=b+1
  }
  pvalue_max=(b-1)/B

  hsphere=TRUE
  b=0
  while (hsphere==TRUE && b<=(B/2)){
    hsphere=(target_est <= Nsample[(B-b)] && Nsample[(b+1)] <= target_est)
    b=b+1

  }
  pvalue_sphere=2*(b-1)/B

  #calculate volumes (lp-Ball)
  vol_max_vec <- rep(NA,length(subsets))
  vol_sphere_vec_upper <- rep(NA,length(subsets))
  vol_sphere_vec_lower <- rep(NA,length(subsets))
  vol_index <- 1
  for (subset in subsets){
    vol_max_vec[vol_index] <- 2^length(subset)*gamma(1/exp+1)^length(subset)/gamma(length(subset)/exp+1)*prod(quant_est[[1]]^(1/exp)*model$sigma_est[subset]/sqrt(n))
    vol_sphere_vec_upper[vol_index] <- 2^length(subset)*gamma(1/exp+1)^length(subset)/gamma(length(subset)/exp+1)*prod(quant_est[[2]]^(1/exp)*model$sigma_est[subset]/sqrt(n))
    vol_sphere_vec_lower[vol_index] <- 2^length(subset)*gamma(1/exp+1)^length(subset)/gamma(length(subset)/exp+1)*prod(quant_est[[3]]^(1/exp)*model$sigma_est[subset]/sqrt(n))
    vol_index <- vol_index + 1
  }

  vol_max <- prod(vol_max_vec)
  vol_sphere <- prod(vol_sphere_vec_upper) - prod(vol_sphere_vec_lower)

  result <- list(estimates = model$estimates,
                 edge_list = model$edge_list,
                 pvalue_max = pvalue_max,
                 pvalue_sphere = pvalue_sphere,
                 hyp_max = hyp_max,
                 hyp_sphere = hyp_sphere,
                 vol_max = vol_max,
                 vol_sphere = vol_sphere,
                 quantiles = quant_est)
  return(result)
}



