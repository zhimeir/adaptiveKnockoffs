#' Adaptive knockoff filter based on Bayeian two group model
#'
#' filter_EM takes the inmportance statistic W as input. The filter determines the order of the statistic by the Bayesian posterior probability of being a null.
#'
#' @param W vector of length p, denoting the imporatence statistics calculated by \code{\link[knockoff]{knockoff.filter}}.
#' @param U p-by-r matrix of side information.
#' @param alpha target FDR level (default is 0.1).
#' @param reveal_prop The proportion of hypotheses revealed at intialization (default is 0.4).
#' @param mute whether \eqn{\hat{fdp}} of each iteration is printed (defalt is TRUE).
#' @param df Degree of freedom of the splines (default is 3).
#' @param R Number of iterations in the EM algorithm
#'
#' @return A list of the following:
#'  \item{nrejs}{The number of rejections for each specified target fdr (alpha) level}
#'  \item{rejs}{Rejsction set fot each specified target fdr (alpha) level}
#'  \item{rej.path}{The order of the hypotheses (used for diagnostics)}
#'  \item{unrevealed.id}{id of the hypotheses that are nor revealed in the end (used for diagnostics)}
#'  \item{tau}{Threshold of each target FDR level (used for diagnostics)}
#'  \item{acc}{The accuracy of classfication at each step (used for diagnostics)}
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
#' #Generate importance statistic using knockoff package
#' W = stat.glmnet_coefdiff(X,Xk,y)
#'
#' #Using filer_EM to obtain the final rejeciton set
#' U = 1:p #Use the location of the hypotheses as the side information
#' result = filter_EM(W,U)
#'
#' @export



filter_EM <- function(W,U,alpha = 0.1,offset = 1,mute = TRUE,df = 3,R=1,s0 = 5e-3,cutoff = NULL){
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
  tau.sel = c()


  # Initializing the output
  rejs = vector("list",length(alpha))
  nrejs =  rep(0,length(alpha))
  ordered_alpha = sort(alpha,decreasing = TRUE)
  rej.path = c()
  index.rank = rep(0,p)

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
      break
    }

  }

  index.rank[unrevealed_id] = length(revealed_id)+1
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,tau = tau.sel,W=W,index.rank = index.rank)
  return(result)
}


## Auxiliary funcitons
calculate.fdphat = function(W,tau,offset = 1){
  p = length(W)
  fdphat = (offset+sum(W<=-tau))/sum(W>=tau)
  return(fdphat)
}


logis <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}

dens  <- function(t,mu,delta){
  val = delta*(t==1/2)+(1-delta)*(t!=1/2)*(t^(1/mu-1)/mu)
  return(val)
}


prob_revealed <- function(pi,mu0,mu1,t,delta0,delta1){
  num = pi*dens(t,mu1,delta1)
  denom = pi*dens(t,mu1,delta1) + (1-pi)*dens(t,mu0,delta0)
  if(denom == 0 ){value =0.5}else{value = num/denom}
  return(value)
}

prob_unrevealed <- function(pi,mu0,mu1,t,delta0,delta1){
  num = pi*dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1)
  denom = num +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if((denom) ==0){value = 0.5}else{value = num/denom}
  return(value)
}


exp_unrevealed_1 <- function(pi,mu0,mu1,t,delta0,delta1){
  y1  = -log(t)
  y2 = -log(1-t)
  num = y1 * pi *dens(t,mu1,delta1) + y2*pi*dens(1-t,mu1,delta1)
  denom =pi *dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1) +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if(denom ==0){value = (y1+y2)/2}else{value = num/denom}
  return(value)
}

exp_unrevealed_0 <- function(pi,mu0,mu1,t,delta0,delta1){
  y1  = -log(t)
  y2 = -log(1-t)
  num = y1 *( 1-pi) *dens(t,mu0,delta0) + y2*(1-pi)*dens(1-t,mu0,delta0)
  denom =pi *dens(t,mu1,delta1) +pi*dens(1-t,mu1,delta1) +(1-pi)*dens(t,mu0,delta0) +(1-pi)*dens(1-t,mu0,delta0)
  if(denom ==0){value = (y1+y2)/2}else{value = num/denom}
  return(value)
}


order_prob <- function(pi,mu0,mu1,t,delta0,delta1){
  h11 = delta1*(t==1/2)+(1-delta1)*(t!=1/2)*(t^{1/mu1-1}/mu1)
  h10 = delta1*(t==1/2)+(1-delta1)*(t!=1/2)*((1-t)^{1/mu1-1}/mu1)
  h01 = delta0*(t==1/2)+(1-delta0)*(t!=1/2)*(t^{1/mu0-1}/mu0)
  h00 = delta0*(t==1/2)+(1-delta0)*(t!=1/2)*((1-t)^{1/mu0-1}/mu0)
  num = pi*h11
  denom = pi*h11+pi*h10 +(1-pi)*(h00+h01)
  if(denom ==0){value =0.5}else{value = num/denom}
  return(value)
}
