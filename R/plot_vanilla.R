#' Diagnostics for adative knockoff results
#'
#' Plot_vanilla plots the ordering of vanilla knockoff as reference.
#'
#' @param W vector of length p, denoting the imporatence statistics calculated by \code{\link[knockoff]{knockoff.filter}}.
#' @param alpha A vector of all target FDR levels.
#' @param nonzero The true signals (default is NULL).
#' @return A plot of the realized ordering by vanilla knockoff.
#'
#' @family plot




plot_vanilla <- function(W,alpha,nonzero = NULL,start_index=1){
  require("ggplot2")
  require("knockoff")

  alpha = sort(alpha,decreasing = TRUE)
  p = length(W)
  tau = c()
  wpath = W[order(abs(W))]
  rej_path = order(abs(W))
  #obtain tau
  for (i in alpha){
    taunew = knockoff.threshold(W,fdr = i,offset = 1)
    tau = c(tau,which(abs(wpath-taunew) == min(abs(wpath-taunew) )))
  }
  if(is.null(nonzero) == 0){
    group = rep(1,p)
    group[nonzero] = 2

    grouppath = group[rej_path]

    df = data.frame(Ordering = start_index:p, W=wpath[start_index:p])
    df$title = "Vanilla Knockoff"
    ggplot(df,aes(Ordering,W))+
      geom_bar(stat = "identity", aes(fill =grouppath[start_index:p]),show.legend = FALSE)+
      geom_vline(xintercept = tau,color = "red", linetype = "dashed",size = 1)+
      theme_bw()+
      theme(axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15),
            axis.text.x=element_text(size = 10),
            axis.text.y=element_text(size = 10),
            panel.border = element_rect(size = 1.5, linetype='solid')
      )

  }else{
    df = data.frame(x = 1:p, W=wpath)
    ggplot(df,aes(x,W))+
      geom_bar(stat = "identity",show.legend = FALSE)+
      geom_vline(xintercept = tau,color = "red")+
      theme_bw()+
      theme(axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15),
            axis.text.x=element_text(size = 10),
            axis.text.y=element_text(size = 10)
      )
  }
}
