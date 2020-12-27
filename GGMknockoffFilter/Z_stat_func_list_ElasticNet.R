### Generate a list of function for constructing Z-statistics. The output is used as an argument for the function 'GKF-Re'
# Here the functions are based on 'Z_coef' and 'Z_max_lambda" with different values of 'alpha' (and 'lambda_quantile_vec')
# The functions in the list should only have arguments 'Nodewise_Y_Xs_Xk_list' and 'num.cores'

##########################################################
### Main function:
Z_stat_func_list_ElasticNet <-function(alpha_vec, lambda_quantile_vec){
  
  list_coef_lambda <- lapply(alpha_vec, function(alpha_used, lambda_quantile_vec) {
    
    force(alpha_used)    # Deal with lazy evaluation
    
    function(Nodewise_Y_Xs_Xk_list, num.cores){
      return(Z_coef(Nodewise_Y_Xs_Xk_list, num.cores, alpha=alpha_used, lambda_quantile_vec))
    }
  }, lambda_quantile_vec)
  
  list_max_lambda <- lapply(alpha_vec, function(alpha_used) {
    
    force(alpha_used)
    
    function(Nodewise_Y_Xs_Xk_list, num.cores){
      
      return(Z_max_lambda(Nodewise_Y_Xs_Xk_list, num.cores, alpha=alpha_used))
    }
  })
  
  list_func <- c(list_coef_lambda, list_max_lambda)
  
  return(list_func)
}
