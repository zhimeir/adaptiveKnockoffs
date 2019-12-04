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


filter_randomForest_getorder <- function(W,z,alpha =0.1,offset=1,reveal_prop = 0.5,mute = TRUE){

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
  alpha = c(alpha,0)
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




filter_EM_getorder <- function(W,U,alpha = 0.1,offset = 1,mute = TRUE,df = 3,R=1,s0 = 5e-3,cutoff = NULL){
  #Check the input format
  if(is.numeric(W)){
    W = as.vector(W)
  }else{
    stop('W is not a numeric vector')
  }

  if(is.numeric(U) ==1){
    U = as.matrix(U)
  }else{
    stop('U is not numeric')
  }

  #if(is.numeric(reveal_prop)==0) stop('reveal_prop should be a numeric.')
  #if(reveal_prop>1) stop('reveal_prop should be a numeric between 0 and 1')
  #if(reveal_prop<0) stop('reveal_prop should be a numeric between 0 and 1')


  #Extract dimensionality
  p = length(W)
  #check if z is in the correct form
  if(dim(U)[1]!=p){
    if(dim(U)[2]==p){
      U = t(U)
    }
    else{
      stop('Please check the dimensionality of the side information!')
    }
  }
  pz = dim(U)[2]
  all_id = 1:p



  # Initializing the output
  rejs = vector("list",length(alpha))
  nrejs =  rep(0,length(alpha))
  alpha = c(alpha,0)
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  index.rank = rep(0,p)
  tau.sel = c()

  # Algorithm initialization
  W_abs = abs(W)
  W_revealed = W_abs


  # Revealing a small amount of signs based on magnitude only
  tau = rep(s0,p)

  revealed_id = which(W_abs<=tau)


  if (length(revealed_id)>0)
  {unrevealed_id =all_id[-revealed_id]
  W_revealed[revealed_id] = W[revealed_id]
  rej.path = c(rej.path,revealed_id)
  index.rank[revealed_id] = 1
  }else{
    unrevealed_id = all_id
  }

  pi = rep(sum(W>0)/p,p)
  delta0 = sum(W==0)/p*(1-mean(pi))
  delta1 = sum(W==0)/p*(mean(pi))
  t = logis(W_revealed)
  mu_0 =  rep(-log(logis(mean(W_abs[W<0]))),p)
  mu_1 = rep(-log(logis(mean(W_abs[W>0]))),p)
  mu_1[W==0] = log(2)
  mu_0[W==0] = log(2)

  H = rep(1e-10,p)
  y0 = -log(t)
  y1 = -log(t)
  count = 0

  for (talpha in 1:length(alpha)){

    fdr = ordered_alpha[talpha]
    for (i in 1:length(unrevealed_id)){
      # EM algorithm
      for (r in 1:R){
        # E step
        # revealed part
        H[revealed_id] = sapply(1:length(revealed_id),function(j) prob_revealed(pi[revealed_id[j]],mu_0[revealed_id[j]],mu_1[revealed_id[j]],t[revealed_id[j]],delta0,delta1))
        y1[revealed_id] = -log(t[revealed_id])
        y0[revealed_id] = -log(t[revealed_id])

        # unrevealed part
        H[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) prob_unrevealed(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))
        H = pmax(H,1e-10)
        H = pmin(H,1-1e-10)
        y1[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) exp_unrevealed_1(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))/H[unrevealed_id]
        y0[unrevealed_id] = sapply(1:length(unrevealed_id),function(j) exp_unrevealed_0(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))/(1-H[unrevealed_id])


        # M step
        if(is.null(cutoff) == TRUE){cutoff=p}
        if(length(unrevealed_id)<cutoff){
          if(dim(U)[2] ==1){
            mdl = gam(log(H/(1-H))~ns(U,df))
            pi = logis(mdl$fitted.values)
          }else{
            mdl = gam(log(H/(1-H))~s(U[,1],U[,2]))
            pi = logis(mdl$fitted.values)
          }
          if(dim(U)[2]==1){mdl =gam(y0[t!=1/2]~ns(U[t!=1/2],df),weights = (1-H[t!=1/2]),family = Gamma(link = "log"))
          }else{
            mdl =gam(y0[t!=1/2]~s(U[t!=1/2,1],U[t!=1/2,2]),weights = (1-H[t!=1/2]),family = Gamma(link = "log"))
          }
          mu_0[t!=1/2] =mdl$fitted.values

          if(dim(U)[2]==1){mdl =gam(y1[t!=1/2]~ns(U[t!=1/2],df),weights = (H[t!=1/2]),family = Gamma(link = "log"))
          }else{
            mdl =gam(y1[t!=1/2]~s(U[t!=1/2,1],U[t!=1/2,2]),weights = (H[t!=1/2]),family = Gamma(link = "log"))
          }
          mu_1[t!=1/2] =mdl$fitted.values
        }else{
          mdl = randomForest(y= H,x=as.matrix(U))
          pi = mdl$predicted
          mdl = randomForest(y = y0[t!=1/2],x= as.matrix(U[t!=1/2,]))
          mu_0[t!=1/2] =mdl$predicted
          mdl = randomForest(y = y1[t!=1/2],x= as.matrix(U[t!=1/2,]))
          mu_1[t!=1/2] =mdl$predicted
        }

        delta0 = sum((1-H)*(t==1/2))/(sum((1-H)*(t==1/2))+sum((1-H)*(t!=1/2)))
        delta1 = sum((H)*(t==1/2))/(sum((H)*(t==1/2))+sum((H)*(t!=1/2)))
      }


      horder = sapply(1:length(unrevealed_id),function(j) order_prob(pi[unrevealed_id[j]],mu_0[unrevealed_id[j]],mu_1[unrevealed_id[j]],t[unrevealed_id[j]],delta0,delta1))
      #ploth = rep(0,p)
      #ploth[unrevealed_id] = sign(W[unrevealed_id])*horder
      #plot(ploth)

      ind.min = which(horder == min(horder))
      if(length(ind.min)==1){
        ind.reveal = ind.min
      }else{
        ind.reveal = ind.min[which.min(W_abs[ind.min])]
      }

      ind.reveal = unrevealed_id[ind.reveal]
      index.rank[ind.reveal] = length(revealed_id)+1
      revealed_id = c(revealed_id,ind.reveal)
      unrevealed_id = all_id[-revealed_id]
      rej.path = c(rej.path,ind.reveal)
      tau[ind.reveal] =  W_abs[ind.reveal]+1
      fdphat = calculate.fdphat(W,tau,offset = offset)
      if(mute == FALSE) print(fdphat)
      if(fdphat<=fdr | fdphat ==Inf){break}
      W_revealed[ind.reveal] = W[ind.reveal]
      t = logis(W_revealed)
    }

    #Check if the estimated FDR is below the target FDR threshold
    if(fdphat<=fdr){
      rej = which(W>=tau)
      rejs[[talpha]] = rej
      nrejs[talpha] = length(rej)
      tau.sel = c(tau.sel, length(revealed_id))
    }else{
      tau.sel = c(tau.sel, length(revealed_id))
      break
    }

  }

  index.rank[unrevealed_id] = length(revealed_id)+1
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,tau = tau.sel,W=W,index.rank = index.rank,tau = tau.sel)
  return(result)
}
