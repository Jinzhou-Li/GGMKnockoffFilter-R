### GGM knockoff filter with fixed hyperparameters (Algprithm 3 in the paper)

# input: 'X': data matrix (n \times p)
#        'q': nominal FDR level
#        'offset': indicate FDR (offset=1) / mFDR (offset=0) control
#        'Z_matrix_function': function used to construct Z-statistics (measure the importance of variable/knockoff to the response)
#        '...': extra arguments of 'Z_matrix_function' besides 'Nodewise_Y_Xs_Xk_list' and 'num.cores'
#        'knockoff_method': choose from ("equi" and "sdo"), method to construct knockoffs
#        'rule': choose from ("AND" and "OR"), rule used to recover the estimated edge set from the estimated neighborhoods
#        'a': choose from (1 and 0.01), parameter used in the optimization problem
#        'W_matrix_BasedOnZ': choose from ("difference" and "max_sign"): method used to construct feature statsitics (i.e. W-statistics) based on Z-statistics
#        'num.cores': number of cores used for the parallel, the default value is 1
# output: 'E_est': estimated edge set

##########################################################

library(knockoff)        # for nodewisely constructing knockoffs
library(parallel)        # parallel when nodewisely construct knockoffs and Z-statistics

##########################################################
### Main function:
GKF <- function(X, q, offset, Z_matrix_function, ..., knockoff_method, rule, a, W_matrix_BasedOnZ, num.cores=1){
  
  # Obtain Y, design matrix X, and X_knockoffs nodewisely
  Nodewise_Y_Xs_Xk_list <- Nodewise_Y_Xs_Xk(X, knockoff_method, num.cores)
  
  # construct Z matrix
  Z_matrix_list <- Z_matrix_function(Nodewise_Y_Xs_Xk_list, num.cores, ...)
  Z_ori_matrix <- Z_matrix_list[[1]]
  Z_knock_matrix <- Z_matrix_list[[2]]
  
  # in case Z_coef_function is used and lambda_quantile is a vector
  if(ncol(Z_ori_matrix)!=nrow(Z_ori_matrix)) {stop("Error: Z_matrix is not a square matrix")}
  
  # construct W matrix
  if(W_matrix_BasedOnZ=="difference"){
    W_matrix <- Z_ori_matrix - Z_knock_matrix
  }
  if(W_matrix_BasedOnZ=="max_sign"){
    W_matrix <- sign(Z_ori_matrix-Z_knock_matrix) * pmax(Z_ori_matrix,Z_knock_matrix)
  }
  
  # compute threshold and return the estimated graph
  E_est <- E_est_givenW_func(W_matrix, q, offset, rule, a)
  
  return(E_est)
}