filter_gam_alpha <- function(W, z, df = 2, alpha = 0.1,offset = 1,
                       reveal_prop = 0.1, mute = TRUE){

  #Check the input format
  if(is.numeric(W)){
    W <- as.vector(W)
  }else{
    stop('W is not a numeric vector.')
  }

  if(is.numeric(z) ==1){
    z <- as.matrix(z)
  }else{
    stop('z is not numeric.')
  }

  if(is.numeric(reveal_prop) == 0) stop('reveal_prop should be a numeric.')
  if(reveal_prop > 1) stop('reveal_prop should be a numeric between 0 and 1.')
  if(reveal_prop < 0) stop('reveal_prop should be a numeric between 0 and 1.')


  #Extract dimensionality
  p <- length(W)
  
  #check if z is in the correct form
  if(dim(z)[1] != p){
    if(dim(z)[2] == p){
      z <- t(z)
    }
    else{
      stop('Please check the dimensionality of the side information!')
    }
  }
  pz <- dim(z)[2]

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
      if(dim(z)[2]==1){
        mdl  =gam(revealed_sign~ns(z,df) + abs(W),family = "binomial")
      }else{
        mdl  =gam(revealed_sign~z + abs(W),family = "binomial")
      }
      fitted.pval = mdl$fitted.values
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

  result = list(rejs = rejs,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,tau = tau.sel,W=W)
  return(result)
}


