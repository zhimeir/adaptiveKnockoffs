#' Adaptive knockoff filter based on GAM (Generalized Additive Model)
#'
#' filter_gam takes the inmportance statistic W as input. 
#' The filter determines the order of the statistic by 
#' the probability of being a null. The probability is 
#' given by GAM.
#'
#' @param W vector of length p, denoting the imporatence statistics calculated by \code{\link[knockoff]{knockoff.filter}}.
#' @param U p-by-r matrix of side information.
#' @param df_list A list of candidates for the degree of freedom of the splines (default is 6:10).
#' @param alpha (a vector of) target FDR levels (default is 0.1).
#' @param reveal_prop The proportion of hypotheses revealed according to |W| (default is 0.1).
#' @param mute whether \eqn{\hat{fdp}} of each iteration is printed (defalt is TRUE).
#'
#' @return A list of the following:
#'  \item{nrejs}{The number of rejections for each specified target fdr (alpha) level}
#'  \item{rejs}{Rejsction set fot each specified target fdr (alpha) level}
#'  \item{rej.path}{The order of the hypotheses (used for diagnostics)}
#'  \item{unrevealed.id}{id of the hypotheses that are nor revealed in the end (used for diagnostics)}
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
#' #Using filer_gam to obtain the final rejeciton set
#' z = 1:p #Use the location of the hypotheses as the side information
#' result = filter_gam(W,z)
#'
#' @export


filter_gam <- function(W, U, alpha = 0.1, offset = 1,
                       df_list = 6:10, reveal_prop = 0.1, 
                       mute = TRUE){

  ## Check the input format
  if(is.numeric(W)){
    W <- as.vector(W)
  }else{
    stop("W should be a numeric vector.")
  }

  if(is.numeric(U)){
    U <- as.matrix(U)
  }else{
    stop("U should be numeric.")
  }

  if(!is.numeric(reveal_prop)) stop("reveal_prop should be a numeric.")
  if(reveal_prop > 1) stop("reveal_prop should be a numeric between 0 and 1.")
  if(reveal_prop < 0) stop("reveal_prop should be a numeric between 0 and 1.")


  ## Extract dimensionality
  p <- length(W)
  
  ## check if U is in the correct form
  if(dim(U)[1] != p){
    if(dim(U)[2] == p){
      U <- t(U)
    }
    else{
      stop("Please check the dimensionality of the side information!")
    }
  }

  ## Prepare the output
  rejs <- list()
  nrejs <- rep(0,length(alpha))
  rej.path <- c()
  
  ## Initilization
  ordered_alpha <- sort(alpha, decreasing = TRUE)
  W_abs <- abs(W)
  W_sign <- (sign(W) + 1) / 2       # code the sign of W_j's as 0, 1, 1/2
  revealed_sign <- rep(1 / 2, p)

  all_id <- 1:p
  revealed_id <- c() 
  unrevealed_id <- all_id
  df_int <- df_list[1]              # set a intial value for df_int in case the model selection fails

  cutoff_val <- quantile(W_abs[W != 0], reveal_prop)
  cutoff <- p - sum(W_abs <= cutoff_val)
  count <- 0

  ## Iteratively reveal hypotheses;
  ## the order is determined by gam()

  for (talpha in 1:length(alpha)){

    fdr <- ordered_alpha[talpha]
    p_remain <- length(unrevealed_id)

    for (i in 1 : p_remain){

      ## If the number of revealed hypothesis is 
      ## less than the cutoff, reveal according 
      ## to |W|

      if(length(unrevealed_id) > cutoff){
        
        ## Reveal the W_j with smallest 
        ## probability of being a positive
        ind.reveal <- which.min(W_abs[unrevealed_id])[1]
      
      }else{
      
        count <- count + 1
      
        if(dim(U)[2] == 1){
        
          ## Select model via BIC every 20 steps
          if(count %% 20 == 1){ 
            ms_res <- select_gam_mdl(df_list, revealed_sign, U, W_abs)
            df_int <- ms_res$opt_df
          }

          mdl <- suppressWarnings(gam(revealed_sign ~ ns(U, df_int) + W_abs, family = binomial()))

        }else{

          mdl <- suppressWarnings(gam(revealed_sign ~ U + W_abs, family = binomial()))

        }
      
        fitted.pval <- mdl$fitted.values[unrevealed_id]

        ## Reveal the W_j with smallest probability of being a positive
        ind.min <- suppressWarnings(which(fitted.pval == min(fitted.pval)))

        if(length(ind.min) == 1){
          ind.reveal <- ind.min
        }else{
          ind.reveal <- ind.min[which.min(W_abs[ind.min])]
        }
      
      }

      ## Update the list of revealed and unrevealed_id
      ind.reveal <- unrevealed_id[ind.reveal]
      revealed_id <- c(revealed_id, ind.reveal)
      unrevealed_id <- all_id[-revealed_id]
      rej.path <- c(rej.path, ind.reveal)
      revealed_sign[ind.reveal] <- W_sign[ind.reveal]
      fdphat <- compute_fdphat(W, revealed_id, offset = offset)

      if(!mute){
        outmsg <- sprintf("The estimated FDP is %.3f.\n", fdphat)
        print(outmsg)
      }

      if(fdphat <= fdr) break

    }

    ## Collect results
    rej <- unrevealed_id[W[unrevealed_id] > 0]
    rejs[[talpha]] <- rej
    nrejs[talpha] <- length(rej)
  }

  result <- list(rejs = rejs, nrejs = nrejs,
                 rej.path = c(rej.path, unrevealed_id),
                 unrevealed_id = unrevealed_id)

  return(result)
}

## Select the optimal degree of freedom via BIC

select_gam_mdl <- function(df_list, signw, U, W_abs){
  
  val <- c()

  for(df in df_list){
   
    gam_mdl <- suppressWarnings(gam(signw  ~ ns(U, df) + W_abs, family = binomial()))
    val <- c(val, BIC(gam_mdl))

  }

  opt_df <- df_list[which.min(val)]

  return(list(opt_df = opt_df))
  
}
