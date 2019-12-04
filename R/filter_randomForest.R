#' Adaptive Knockoff Filter With Random Forest
#'
#' filter_glm returns a set of rejections with FDR controlled at custom target
#'
#' @param W vector of length p, denoting the imporatence statistics calculated by \code{\link[knockoff]{knockoff.filter}}.
#' @param z p-by-r matrix of side information.
#' @param alpha target FDR level (default is 0.1).
#' @param offset either 0 or 1 (default: 1). The offset used to compute the rejection threshold on the statistics. For details, see \code{\link[knockoff]{knockoff.threshold}}.
#' @param reveal_prop The proportion of hypotheses revealed at intialization (default is 0.5).
#' @param mute whether \eqn{\hat{fdp}} of each iteration is printed (defalt is TRUE).
#'
#' @return A list of the following:
#'  \item{nrejs}{The number of rejections for each specified target fdr (alpha) level}.
#'  \item{rejs}{Rejsction set fot each specified target fdr (alpha) level}.
#'  \item{rej.path}{The order of the hypotheses (used for diagnostics)}.
#'  \item{unrevealed.id}{id of the hypotheses that are nor revealed in the end (used for diagnostics)}.
#'  \item{tau}{Threshold of each target FDR level (used for diagnostics)}.
#'  \item{acc}{The accuracy of classfication at each step (used for diagnostics)}.
#'
#'
#'
#' @family filter
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
#' #Gnerate importance statistic using knockoff package
#' W = stat.glmnet_coefdiff(X,Xk,y)
#'
#' #Using filer_gam to obtain the final rejeciton set
#' z = 1:p #Use the location of the hypotheses as the side information
#' result = filter_glm(W,z)
#'

#' @export


filter_randomForest <- function(W,z,alpha =0.1,offset=1,reveal_prop = 0.5,mute = TRUE){

  #Check the input format
  if(is.numeric(W)){
    W = as.vector(W)
  }else{
    stop('W is not a numeric vector')
  }

  if(is.numeric(z) ==1){
    z = as.matrix(z)
  }else{
    stop('z is not numeric')
  }

  if(is.numeric(reveal_prop)==0) stop('reveal_prop should be a numeric.')
  if(reveal_prop>1) stop('reveal_prop should be a numeric between 0 and 1')
  if(reveal_prop<0) stop('reveal_prop should be a numeric between 0 and 1')


  #Extract dimensionality
  p = length(W)
  #check if z is in the correct form
  if(dim(z)[1]!=p){
    if(dim(z)[2]==p){
      z = t(z)
    }
    else{
      stop('Please check the dimensionality of the side information!')
    }
  }
  pz = dim(z)[2]

  #Initilization
  tau = rep(0,p)
  rejs = list()
  nrejs = rep(0,length(alpha))
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  W_abs = abs(W)
  W_sign = as.numeric(W>0)
  revealed_sign = rep(1,p)
  all_id = 1:p
  tau.sel = c()
  acc = c()


  #Reveal a small proportion of W
  revealed_id = which(W_abs<=quantile(W_abs,reveal_prop))

  #Update the revealed information
  revealed_sign[revealed_id] = W_sign[revealed_id]
  unrevealed_id =all_id[-revealed_id]
  tau[revealed_id] = W_abs[revealed_id]+1
  rej.path = c(rej.path,revealed_id)

  #Iteratively reveal hypotheses; the order determined by glmnet
  for (talpha in 1:length(alpha)){

    fdr = ordered_alpha[talpha]

    for (i in 1:length(unrevealed_id)){

      mdl = randomForest(y = as.factor(revealed_sign),x = cbind(W_abs,z),norm.votes = TRUE,ntree = 1000)
      fitted.pval = mdl$votes[,ncol(mdl$votes)]
      fitted.pval = fitted.pval[unrevealed_id]
      predicted.sign = fitted.pval>0.5
      acc = c(acc, sum(predicted.sign == W_sign[unrevealed_id])/length(unrevealed_id))

      #Reveal the W_j with smallest probability of being a positive
      ind.min = which(fitted.pval == min(fitted.pval))
      if(length(ind.min)==1){
        ind.reveal = ind.min
      }else{
        ind.reveal = ind.min[which.min(W_abs[ind.min])]
      }

      ind.reveal = unrevealed_id[ind.reveal]
      revealed_id = c(revealed_id,ind.reveal)
      rej.path = c(rej.path,ind.reveal)
      unrevealed_id = all_id[-revealed_id]
      revealed_sign[ind.reveal] = W_sign[ind.reveal]
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W,tau,offset = offset)
      if (mute == "FALSE"){print(fdphat)}
      if(fdphat<=fdr){break}
    }
    rej = which(W>=tau)
    rejs[[talpha]] = rej
    nrejs[talpha] = length(rej)
    tau.sel = c(tau.sel, length(revealed_id))
  }

  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,tau = tau.sel,acc = acc,W=W,alpha = alpha)
  return(result)
}


