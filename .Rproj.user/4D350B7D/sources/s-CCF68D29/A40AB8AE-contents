#' Create the confidence intervals for a GGMtest object
#' @export

confint.GGMtest <- function(object){
  checkmate::check_choice(object$additional_parameters$s,1)
  CI <- data.frame(object$estimates - object$additional_parameters$sigma_est*object$quantile[1]/sqrt(object$additional_parameters$n),
          object$estimates + object$additional_parameters$sigma_est*object$quantile[1]/sqrt(object$additional_parameters$n))
  colnames(CI) <- c(sprintf(" %s%%",object$additional_parameters$alpha/2*100),sprintf(" %s%%",(1-object$additional_parameters$alpha/2)*100))
  rownames(CI) <- apply(object$edge_list, 1, function(x) sprintf("(%i,%i)", x[1], x[2]))
  return(CI)
}




