### Generate a list of function for constructing Z-statistics. The output is used as an argument for the function 'GKF-Re'
# Here the functions are based on 'Z_coef' and 'Z_max_lambda" with different values of 'alpha' (and 'lambda_quantile_vec')
# The functions in the list should only have arguments 'Nodewise_Y_Xs_Xk_list' and 'num.cores'

##########################################################
### Main function:
Z_stat_func_list_ElasticNet <-function(alpha_vec, lambda_quantile_vec){
  
  func.list <- list()
  
  for (i in 1:length(alpha_vec)) {
    func.list[[i]] <- function(Nodewise_Y_Xs_Xk_list, num.cores){
      return(Z_coef(Nodewise_Y_Xs_Xk_list, num.cores, alpha=alpha_vec[i], lambda_quantile_vec))
    }
  }
  
  for (i in 1:length(alpha_vec)) {
    func.list[[length(alpha_vec)+i]] <- function(Nodewise_Y_Xs_Xk_list, num.cores){
      return(Z_max_lambda(Nodewise_Y_Xs_Xk_list, num.cores, alpha=alpha_vec[i]))
    }
  }
  
  return(func.list)
}
