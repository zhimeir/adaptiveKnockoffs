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



plot_rank <- function(mdl,nonzero = NULL){
  #require ggplot2 package
  require("ggplot2")
  
  index.rank = mdl$index.rank
  W = mdl$W
  p = length(W)
  
  if(is.null(nonzero) == 0){
    nzid = which(W!=0)
    group = rep("gray",p)
    group[nonzero] = "royalblue"
    group = group[nzid]
    df = data.frame(x = 1:length(nzid), W=W[nzid])
    pp = ggplot(df,aes(x,W))+
      geom_bar(stat = "identity",show.legend = FALSE,aes(fill = group,alpha=0.7))+
     geom_point(aes(x=x,y=index.rank[nzid]*max(W)/max(index.rank)),color = "gray50",size=2)+
      #geom_smooth(aes(x=x,y=index.rank[nzid]*max(W)/max(index.rank)),color = "black",size=0.5,se = FALSE)+
      theme_bw()
      print(pp)
  }else{
    df = data.frame(x = start_index:p, W=wpath[start_index:p])
    pp = ggplot(df,aes(x,W))+
      geom_bar(stat = "identity",show.legend = FALSE)+
      geom_point(aes(x=1:p,y=index.rank*max(W)/max(index.rank)),color = "red")+
      theme_bw()
  }
  
  return(pp)
}
