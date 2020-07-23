### Z-statistics constructed based on the estimated coefficient using elastic net

library(glmnet)
# library(parallel)

# Main function:
Z_coef <- function(Nodewise_Y_Xs_Xk_list, num.cores, alpha, lambda_quantile_vec){
  
  Z_list <- mcmapply(Z_coef_OneNode, 1:ncol(Nodewise_Y_Xs_Xk_list), 
                   MoreArgs=list(Nodewise_Y_Xs_Xk_list, alpha, lambda_quantile_vec), mc.cores=num.cores)
  
  # when 'length(lambda_quantile_vec)>1', there are many groups of hyperparameters, so return many Z_matrix.
  num_Z_matrix <- length(lambda_quantile_vec)
  Z_ori_matrix_array <- Z_knock_matrix_array <- array(0, dim=c(ncol(Nodewise_Y_Xs_Xk_list), ncol(Nodewise_Y_Xs_Xk_list), num_Z_matrix))
  
  for (i in 1:num_Z_matrix) {
    Z_ori_matrix_array[,,i] <- do.call(rbind, lapply(Z_list[1,], `[`,i,))
    Z_knock_matrix_array[,,i] <- do.call(rbind, lapply(Z_list[2,], `[`,i,))
  }
  
  return(list(Z_ori_matrix_array, Z_knock_matrix_array))
}

### sub-function
Z_coef_OneNode <- function(node_index, Nodewise_Y_Xs_Xk_list, alpha, lambda_quantile_vec){
  
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
  
  # generate Z-statistics for orginal variables and knockoffs
  fit <- glmnet(cbind(Xs.swap, Xk.swap), Y, alpha=alpha, intercept=FALSE)
  lambda_vector <- fit$lambda
  lambda_select <- lambda_vector[ceiling(lambda_quantile_vec*length(lambda_vector))]
  # coefficient matrix, each column corresponds to a different 'lambda_quantile' value, row is the coef of each variable
  coef_matrix <- as.matrix(coef(fit, s=lambda_select)[-1,])
  
  Z_ori.swap <- t(abs(coef_matrix[1:dim.variable, ]))
  Z_knock.swap <- t(abs(coef_matrix[(dim.variable+1):(2*dim.variable), ]))
  
  # Correct for swapping of columns of Xs and Xk
  swap.M2 = matrix(swap, nrow=nrow(Z_ori.swap), ncol=length(swap), byrow=TRUE)
  Z_ori <- Z_ori.swap * (1-swap.M2) + Z_knock.swap * swap.M2
  Z_knock <- Z_ori.swap * swap.M2 + Z_knock.swap * (1-swap.M2)
  
  Z_ori_add0 <- Z_knock_add0 <- matrix(0, nrow=nrow(Z_ori), ncol=ncol(Xs)+1)
  Z_ori_add0[,-node_index] <- Z_ori
  Z_knock_add0[,-node_index] <- Z_knock
  
  return(list(Z_ori_add0, Z_knock_add0))
}
