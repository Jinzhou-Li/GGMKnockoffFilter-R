### Return a list of response vector 'Y', standardized (column length 1) design matrix 'Xs' and knockoff matrix 'Xk'
# Correpsond linear models obtained by treating each node as response

# library(parallel)

Nodewise_Y_Xs_Xk <- function(X, knockoff_method, num.cores){
  
  Nodewise_Y_Xs_Xk_list <- mcmapply(OneNode_func, 1:ncol(X), MoreArgs=list(X, knockoff_method), mc.cores=num.cores)
  return(Nodewise_Y_Xs_Xk_list)
}

# Sub-function: Treat one node as response 
# Return response vector, standardized design matrix and knockoff matrix
OneNode_func <- function(node_index, X, knockoff_method){
  
  X_design <- X[, -node_index]
  Y <- X[,node_index]
  
  X_Xk <- create.fixed(X_design, method=knockoff_method)   # function from "knockoffs" package
  Xs <- X_Xk$X
  Xk <- X_Xk$Xk
  
  return(list(Y, Xs, Xk))
}
