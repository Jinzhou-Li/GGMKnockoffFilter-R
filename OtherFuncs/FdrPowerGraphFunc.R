####################################################################################################################
############ Calculate fdp and power based on the estimated edge set(s) E and true precision matrix Omega #########

#######################
##### For one estimated E
Fdp_Power_Graph_func <- function(E_est, Omega){
  
  if(isSymmetric(E_est)==FALSE) {print("E_est is not symmetric!")}
  if(sum(abs(diag(E_est)))!=0) {print("Diagonals of E_est are not 0!")}
  
  ### TRUE adjacent matrix 
  diag(Omega) <- 0
  adj.true <- (Omega != 0)+0
  num.edge <- sum(adj.true > 0) /2    # number of true edges
  
  ### calculate fdp and power
  num.dis <- sum(E_est > 0) /2                  # number of discoveries
  num.fd <- sum((adj.true - E_est) < 0)/2       # number of false discoveries
  num.td <- sum( 2*adj.true - E_est == 1 )/2    # number of true discoveries
  
  fdp <- num.fd / max(num.dis,1)
  if(num.edge==0) {power <- 0} else {power <- num.td/num.edge}                    
  
  return( c(fdp, power) )
}

#######################
##### For a list of estimated E
Fdp_Power_GraphList_func <- function(E_est_list, Omega){
  
  num_E <- length(E_est_list)
  
  fdp_power_index_func <- function(index,E_est_list, Omega){
    E_est <- E_est_list[[index]]
    return(Fdp_Power_Graph_func(E_est, Omega))
  }
  
  fdr_power_all_algo <- t(mapply(fdp_power_index_func, 1:num_E, MoreArgs=list(E_est_list, Omega)))
  return(fdr_power_all_algo)
}
