#' Create adjacency matrix for GGMtest object
#'
#' @param object GGMtest object
#'
#' @return An adjacency matrix with TRUE for all rejected edges. All other edges are set to FALSE.
#' @seealso \code{\link{GGMtest}}
#'
#' @export


adj_GGMtest <- function(object){
  checkmate::checkClass(object, "GGMtest")
  CI.data <- stats::confint(object)
  CI.data[,"Hypothesis"] <- 0>CI.data[,1]&CI.data[,2]>0
  rejected_hyp <- object$edge_list[!CI.data[,"Hypothesis"],]

  adj_mat <- matrix(F,nrow = object$additional_parameters$p, ncol = object$additional_parameters$p)
  invisible(apply(rejected_hyp,1,function(x) adj_mat[x[1],x[2]] <<- T))
  invisible(apply(rejected_hyp,1,function(x) adj_mat[x[2],x[1]] <<- T))
  return(adj_mat)
}


