#### Generate 'Z_stat_function_list' based on elastic net (coef and max_lambda) for 'GKF-Re'

# Z-function in this list should only has argument 'Nodewise_Y_Xs_Xk_list' and 'num.cores' left

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