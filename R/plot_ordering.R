#' Diagnostics for adative knockoff results
#'
#' Plot_ordering takes the adaptive knockoff result as input and plot the ordering.
#'
#' @param mdl the result given by adaptive knockoff filters.
#' @param nonzero The true signals (default is NULL).
#' @return A plot of the realized ordering of hypotheses.
#'
#' @family plot
#'
#' @examples
#' #Generating data
#' p=100;n=100;k=40;
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = 1:k
#' beta = 5*(1:p%in%nonzero)*sign(rnorm(p))/ sqrt(n)
#' y = X%*%beta + rnorm(n,1)
#'
#' #Generate knockoff copy
#' Xk = create.gaussian(X,mu,Sigma)
#'
#' #Generate importance statistic using knockoff package
#' W = stat.glmnet_coefdiff(X,Xk,y)
#'
#' #Using filer_EM to obtain the final rejeciton set
#' U = 1:p #Use the location of the hypotheses as the side information
#' result = filter_randomForest(W,U)
#' plot_ordering(result,nonzero = nonzero)
#' @export



plot_ordering <- function(mdl,nonzero = NULL,start_index = 1,alpha=NULL){
  #require ggplot2 package
  require("ggplot2")

  rej_path = mdl$rej.path
  W = mdl$W
  p = length(W)

  wpath = W[rej_path]
  if(is.null(nonzero) == 0){
    group = rep(1,p)
    group[nonzero] = 2

    grouppath = group[rej_path]

    df = data.frame(Ordering = start_index:p, W=wpath[start_index:p])
    pp = ggplot(df,aes(Ordering,W))+
      geom_bar(stat = "identity", aes(fill =grouppath[start_index:p]),show.legend = FALSE)+
      geom_vline(xintercept = mdl$tau[-length(mdl$tau)],color = "red", linetype = "dashed",size = 1)+
      theme_bw()+
      theme(axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15),
            axis.text.x=element_text(size = 10),
            axis.text.y=element_text(size = 10),
            panel.border = element_rect(size = 1.5, linetype='solid')
              )
  }else{
    df = data.frame(Ordering = start_index:p, W=wpath[start_index:p])
    pp = ggplot(Ordering,aes(x,W))+
      geom_bar(stat = "identity",show.legend = FALSE)+
      geom_vline(xintercept = mdl$tau,color = "red")+
      theme_bw()+
      theme(axis.title.x=element_text(size = 15),
            axis.title.y=element_text(size = 15),
            axis.text.x=element_text(size = 10),
            axis.text.y=element_text(size = 10),
            panel.border = element_rect(size = 1.5, linetype='solid')
      )
  }

return(pp)
}
