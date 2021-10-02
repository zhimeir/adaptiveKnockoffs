#' @export
filter_bin_EM <- function(W, U, alpha = 0.1, offset = 1, df = 3, df_list = 1:10,
                      mute = TRUE, reveal_prop = 0.1, R=1, 
                      cutoff = NULL, tol = 1e-4){
  
  ## Check the input format
  if(is.numeric(W)){
    W <- as.vector(W)
  }else{
    stop('W is not a numeric vector')
  }

  if(is.numeric(U) == 1){
    U <- as.matrix(U)
  }else{
    stop('U is not numeric')
  }

  if(is.numeric(reveal_prop) == 0) stop('reveal_prop should be a numeric.')
  if(reveal_prop > 1 | reveal_prop < 0) stop('reveal_prop should be a numeric between 0 and 1')


  ## Extract dimensionality
  p <- length(W)

  ## check if z is in the correct form
  if(dim(U)[1]!=p){
    if(dim(U)[2]==p){
      U <- t(U)
    }
    else{
      stop('Please check the dimensionality of the side information!')
    }
  }
  pz <- dim(U)[2]
  all_id <- 1:p
  tau.sel <- c()


  # Initializing the output
  rejs <- vector("list",length(alpha))
  nrejs <- rep(0,length(alpha))
  ordered_alpha <- sort(alpha,decreasing = TRUE)
  rej.path <- c()
  index.rank <- rep(0,p)

  # Algorithm initialization
  W_abs <- abs(W)
  W_revealed <- W_abs


  # Revealing a small amount of signs based on magnitude only
  s0 <- quantile(W_abs[W_abs > 0], reveal_prop)
  revealed_id <- which(W_abs <= s0)

  if (length(revealed_id)>0){
    unrevealed_id <- all_id[-revealed_id]
    W_revealed[revealed_id] <- W[revealed_id]
    rej.path <- c(rej.path, revealed_id)
    index.rank[revealed_id] <- 1
  }else{
    unrevealed_id <- all_id
  }

 

  ## Initialization
  ## Probability of nulls: nu
  mdl <- gam(1 - (all_id %in% revealed_id) ~ ns(U, df),  family = binomial())
  nu <- mdl$fitted.values
  
  ## Signs: expected value of sign(W)
  S <- rep(0, p)
  S[revealed_id] <- 1 - (1 + sign(W_revealed[revealed_id])) / 2
  
  ## Indicator: expected value of being a null
  H <- nu

  ## eta
  eta <- rep(1 / 2, p)
  mdl <- gam(S[W != 0] ~ ns(U[W != 0], df),  family = binomial())
  eta[W!=0] <- mdl$fitted.values
  count <- 0
  for (talpha in 1:length(alpha)){

    fdr <- ordered_alpha[talpha]

    for (i in 1:length(unrevealed_id)){
      count <- count + 1 
      ## EM algorithm
      for (r in 1:R){ 
        signw <-  1 - (sign(W) + 1) / 2
        ## E step
        ## revealed part
        H[revealed_id] <- sapply(1:length(revealed_id), function(j) 
                                  H_rev(nu[revealed_id[j]], eta[revealed_id[j]], signw[revealed_id[j]]))

        S[revealed_id] <- signw[revealed_id]

        ## unrevealed part
        H[unrevealed_id] <- nu[unrevealed_id]

        S[unrevealed_id] <-sapply(1:length(unrevealed_id), function(j) 
                                  S_unrev(nu[unrevealed_id[j]], eta[unrevealed_id[j]], W_revealed[unrevealed_id[j]]))


        # M step
        if(count %% 20 == 1){
          ## Model selection
          res <- select_mdl(df_list, U, W, H, S)
          ##         res <- cv_select_mdl(df_list, U, W, H, S)
          nu_df <- res$opt_nu_df
          eta_df <- res$opt_eta_df

          if(dim(U)[2] ==1){
            mdl <- gam(H  ~ ns(U, nu_df), family = quasibinomial())
            if(max(abs(mdl$fitted.value - nu)) > tol){
              nu <- mdl$fitted.values
            }
          }else{
            mdl <- gam(log(H / (1-H)) ~ s(U[,1],U[,2]))
            nu <- logis(mdl$fitted.values)
          }
          
          if(dim(U)[2]==1){
            ##             mdl <- gam(S[W != 0] ~ ns(U[W != 0], eta_df) + W_abs[W != 0], family = binomial())
            mdl <- gam(S ~ ns(U, eta_df) + W_abs, family = quasibinomial())
          }else{
            mdl <- gam(y0[t!=1/2]~s(U[t!=1/2,1],U[t!=1/2,2]),weights = (1-H[t!=1/2]),family = Gamma(link = "log"))
          }
          if(max(abs(mdl$fitted - eta[W!=0])) > tol){
            eta[W != 0] <- mdl$fitted.values[W!=0]
          }
        }
      }


      horder <- sapply(1:length(unrevealed_id),function(j) 
                      posterior(nu[unrevealed_id[j]],eta[unrevealed_id[j]],W_revealed[unrevealed_id[j]])
                      )
      

      ind.min <- which(horder == min(horder))
      if(length(ind.min) == 1){
        ind.reveal <- ind.min
      }else{
        ind.reveal <- ind.min[which.min(W_abs[ind.min])]
      }

      ind.reveal <- unrevealed_id[ind.reveal]
      index.rank[ind.reveal] <- length(revealed_id)+1
      revealed_id <- c(revealed_id,ind.reveal)
      unrevealed_id <- all_id[-revealed_id]
      W_revealed[ind.reveal] = W[ind.reveal]
      rej.path <- c(rej.path, ind.reveal)
      fdphat <- compute_fdphat(W, revealed_id, offset = offset)
      
      if(mute == FALSE) print(fdphat)
      if(fdphat<=fdr | fdphat ==Inf){break}
    }

    #Check if the estimated FDR is below the target FDR threshold
    if(fdphat <= fdr){
      rej <- unrevealed_id[which(W[unrevealed_id] > 0)]
      rejs[[talpha]] <- rej
      nrejs[talpha] <- length(rej)
    }else{
      break
    }

  }

  index.rank[unrevealed_id] = length(revealed_id)+1
  result = list(rejs = rejs,fdphat = fdphat,nrejs = nrejs,rej.path = c(rej.path,unrevealed_id),unrevealed_id = unrevealed_id,tau = tau.sel,W=W,index.rank = index.rank)
  return(result)
}


## Auxiliary funcitons
compute_fdphat = function(W, revealed_id, offset = 1){
  p <- length(W)
  fdphat <- (offset + sum(W[-revealed_id]< 0)) / max(sum(W[-revealed_id] > 0), 1)
  return(fdphat)
}


logis <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}

H_rev <- function(nu_id, eta_id, s_id){
 
 if(s_id == 1){
  hexp <-  nu_id * eta_id / (nu_id * eta_id + (1 - nu_id) / 2)
 }

 if(s_id == 0){
   hexp <- nu_id * (1 - eta_id) / (nu_id * (1 - eta_id) + (1 - nu_id) / 2)
 }

 if(s_id == 1/2){
  hexp <- nu_id
 }

 return(hexp)

}

S_unrev <- function(nu_id, eta_id, w_id){

  if(w_id != 0){
    sexp <- nu_id * eta_id + (1 - nu_id) / 2
  }else{
    sexp <- 1 / 2
  }

  return(sexp)
}

posterior <- function(nu_id, eta_id, w_id){
  
  if(w_id == 0){
    prob <- 0
  }else{
    prob <- nu_id * (1 - eta_id)
  }

  return(prob)
}

select_mdl <- function(df_list, U, W, H, S){
  
  nu_val <- c()
  eta_val <- c()

  for(df in df_list){
    nu_mdl <- gam(H  ~ ns(U, df), family = binomial())
    eta_mdl <- gam(S ~ ns(U, df) + abs(W), family = binomial()) 
    nu_val <- c(nu_val, BIC(nu_mdl))
    eta_val <- c(eta_val, BIC(eta_mdl))
  }
  opt_nu_df <- df_list[which.min(nu_val)]
  opt_eta_df <- df_list[which.min(eta_val)]

  return(list(opt_nu_df = opt_nu_df, opt_eta_df = opt_eta_df))
}

cv_select_mdl <- function(df_list, U, W, H, S, nfold = 2){
  
  nu_val <- c()
  eta_val <- c()
  ntest <- floor(p / nfold)
  all_id <- 1:p

  for(fold in nfold){
    test_id <- sample(1:p, ntest)
    train_id <- all_id[-test_id]

    for(df in df_list){
      df_train_nu <- data.frame(H = H[train_id], U = U[train_id], W = abs(W[train_id]))
      df_train_eta <- data.frame(S = S[(all_id %in% train_id)],
                                 U = U[(all_id %in% train_id)], 
                                 W = abs(W[(all_id %in% train_id)]))
      nu_mdl <- gam(H ~ ns(U, df) + W, family = binomial(), data = df_train_nu)
      eta_mdl <- gam(S ~ ns(U, df) + W, family = binomial(), data = df_train_eta) 
  
      df_test <- data.frame(U = U[test_id], W = abs(W[test_id]))
      nu_train <- predict(nu_mdl, df_test)
      eta_train <- predict(eta_mdl, df_test)

      lik_nu <- H[test_id] * log(nu_test) + (1 - H[test_id]) * log(1 - nu_test)
      lik_eta <- H[test_id] * (S[test_id] * log(eta_test) + (1 - S[test_id]) * log(1 - eta_test))
      nu_val <- c(nu_val, lik_nu)
      eta_val <- c(eta_val, lik_eta)
    }
  }
  opt_nu_df <- df_list[which.max(nu_val)]
  opt_eta_df <- df_list[which.max(eta_val)]

  return(list(opt_nu_df = opt_nu_df, opt_eta_df = opt_eta_df))
}
