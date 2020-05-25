#' Plot confidence intervals for GGMtest object
#' @export

plot.GGMtest <- function(object, edges = object$edge_list[1,]){
  if (object$additional_parameters$s == 1){
    edge_sublist <- apply(edges, 1, function(x) sprintf("(%i,%i)", x[1], x[2]))
    CI.data <- na.omit(confint(object)[edge_sublist,])
    names(CI.data) <- c("lower", "upper")
    CI.data[,"point_estimate"] <- apply(CI.data,1,mean)
    CI.data[,"Edge"] <- rownames(CI.data)
    CI.data[,"Hypothesis"] <- 0>CI.data[,1]&CI.data[,2]>0

    Intervals <- ggplot2::ggplot(CI.data,ggplot2::aes(x=Edge, y=point_estimate)) +
      ggplot2::geom_hline(yintercept = 0, colour = c("red")) +
      ggplot2::geom_errorbar(data=CI.data[!CI.data[,"Hypothesis"],],ggplot2::aes(ymin=lower, ymax=upper),size=1,  colour="red",width = .5) +
      ggplot2::geom_errorbar(data=CI.data[CI.data[,"Hypothesis"],],ggplot2::aes(ymin=lower, ymax=upper),size=1,  colour="blue",width = .5) +
      ggplot2::geom_point(size=2, col = "black") +
      ggplot2::labs(x ="Edge", y = "Estimator")
  } else if (object$additional_parameters$s == 2 & object$additional_parameters$exponent == 2){
    # obtain the indices of the desired edges
    edges_vec <- c(apply(edges, 1, function(x) sprintf("(%i,%i)", x[1], x[2])),apply(edges, 1, function(x) sprintf("(%i,%i)", rev(x)[1], rev(x)[2])))
    all_edges_vec <- apply(object$edge_list, 1, function(x) sprintf("(%i,%i)", x[1], x[2]))
    edges_ind <- na.omit(match(edges_vec,all_edges_vec))

    # obtain the random assigned pairs
    random_pairs <- matrix(c(object$additional_parameters$random_order,c(object$additional_parameters$random_order[-1],object$additional_parameters$random_order[1])),ncol = 2,byrow = F)

    # obtain the selected pairs
    selected_pairs <- random_pairs[sort(c(match(edges_ind,random_pairs[,1]),match(edges_ind,random_pairs[,2]))),]

    df <- data.frame("Var_1" = object$estimates[selected_pairs[,1]],
                     "Var_2" = object$estimates[selected_pairs[,2]],
                     "sigma_1" = object$additional_parameters$sigma_est[selected_pairs[,1]]/sqrt(object$additional_parameters$n)*sqrt(object$quantile[1]),
                     "sigma_2" = object$additional_parameters$sigma_est[selected_pairs[,2]]/sqrt(object$additional_parameters$n)*sqrt(object$quantile[1]),
                     "edge_pair" = apply(selected_pairs, 1, function(x) sprintf("(%i,%i) - (%i,%i)", object$edge_list[x[1],][1],object$edge_list[x[1],][2],object$edge_list[x[2],][1],object$edge_list[x[2],][2])),
                     "hypothsis" = cbind(object$additional_parameters$beta_0[selected_pairs[,1]],object$additional_parameters$beta_0[selected_pairs[,2]]))

    Intervals <- ggplot2::ggplot(df) +
      ggforce::geom_ellipse(ggplot2::aes(x0=Var_1,y0=Var_2,a=sigma_1,b=sigma_2,angle = 0,fill = "Confidence Interval"), alpha = 0.2) +
      ggplot2::geom_point(ggplot2::aes(x=Var_1,y=Var_2,color = "Estimated"),alpha = 10) +
      ggplot2::geom_point(ggplot2::aes(x=hypothsis.1,y=hypothsis.2,color = "Hypothesis"),alpha = 1) +
      ggplot2::labs(x="First Edge",y="Second Edge", fill="",color="Precision Matrix") +
      ggplot2::scale_fill_manual(values = c("blue")) +
      ggplot2::scale_color_manual(values = c("blue","red")) +
      ggplot2::facet_wrap(. ~ edge_pair) +
      ggplot2::theme(legend.position="bottom")
  } else{
    stop("Not defined for the specified confidence interval.")
  }
  return(Intervals)
}
