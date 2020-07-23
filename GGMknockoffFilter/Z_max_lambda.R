### Z-statistics constructed based on the maximal lambda enter the path using elastic net

library(glmnet)

##########################################################
### Main function:
Z_max_lambda <- function(Nodewise_Y_Xs_Xk_list, num.cores, alpha){
  Z_list <- mcmapply(Z_max_lambda_OneNode, 1:ncol(Nodewise_Y_Xs_Xk_list), 
                     MoreArgs=list(Nodewise_Y_Xs_Xk_list, alpha), mc.cores=num.cores)
  
  Z_ori_matrix <- do.call(rbind, Z_list[1,])
  Z_knock_matrix <- do.call(rbind, Z_list[2,])
  
  return(list(Z_ori_matrix,Z_knock_matrix))
}

### sub-function 1:
Z_max_lambda_OneNode <- function(node_index, Nodewise_Y_Xs_Xk_list, alpha){
  
  Y <- (Nodewise_Y_Xs_Xk_list[1,])[[node_index]]
  Xs <- (Nodewise_Y_Xs_Xk_list[2,])[[node_index]]
  Xk <- (Nodewise_Y_Xs_Xk_list[3,])[[node_index]]
  
  dim.variable <- ncol(Xs)
  
  # Randomly swap columns of X_n and X_k 
  # (Matteo: to avoid a computational issue of Lasso and compute secure test statistics. 
  # Same trick used in the "knockoffs" package)
  swap = rbinom(ncol(Xs), 1, 0.5)
  swap.M = matrix(swap, nrow=nrow(Xs), ncol=length(swap), byrow=TRUE)
  Xs.swap  = Xs * (1-swap.M) + Xk * swap.M
  Xk.swap = Xs * swap.M + Xk * (1-swap.M)
  
  Z_ori_knock.swap <- lasso_max_lambda_glmnet(cbind(Xs.swap, Xk.swap), Y, alpha=alpha)   # function from "knockoffs" package
  Z_ori.swap <- Z_ori_knock.swap[1:dim.variable]
  Z_knock.swap <- Z_ori_knock.swap[(dim.variable+1):(2*dim.variable)]
  
  # Correct for swapping of columns of Xs and Xk
  Z_ori <- Z_ori.swap * (1-swap) + Z_knock.swap * swap
  Z_knock <- Z_ori.swap * swap + Z_knock.swap * (1-swap)
  
  Z_ori_add0 <- rep(0, ncol(Xs)+1)
  Z_knock_add0 <- rep(0, ncol(Xs)+1)
  
  Z_ori_add0[-node_index] <- Z_ori
  Z_knock_add0[-node_index] <- Z_knock
  
  return(list(Z_ori_add0, Z_knock_add0))
}

### sub-function 2: used to calculate the maximal-lambda enter the lasso path. 
### This function is taken from R-package "knockoff"

lasso_max_lambda_glmnet <- function(X, y, nlambda=500, intercept=T, standardize=T, ...) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  
  # Standardize the variables
  if( standardize ){
    X = scale(X)
  }
  
  n = nrow(X); p = ncol(X)
  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family
  
  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  
  fit <- glmnet::glmnet(X, y, lambda=lambda, intercept=intercept, 
                        standardize=F, standardize.response=F, ...)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  indices <- apply(fit$beta, 1, first_nonzero)
  names(indices) <- NULL
  ifelse(is.na(indices), 0, fit$lambda[indices] * n)
}