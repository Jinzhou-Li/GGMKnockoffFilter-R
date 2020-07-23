### Obtain the related response vector, standardized (column length 1) design matrix and the knockoff matrix of the linear model by treating each node as response

# input: 'X': data matrix (n \times p)
#        'knockoff_method': choose from ("equi" and "sdo"), method to construct knockoffs
#        'num.cores': number of cores used for the parallel
# output: 'Nodewise_Y_Xs_Xk_list': a list of response vector 'Y', standardized design matrix 'Xs' and knockoff matrix 'Xk'

##########################################################

# library(knockoff)
# library(parallel)      

##########################################################
### Main function:
Nodewise_Y_Xs_Xk <- function(X, knockoff_method, num.cores){
  
  Nodewise_Y_Xs_Xk_list <- mcmapply(OneNode_func, 1:ncol(X), MoreArgs=list(X, knockoff_method), mc.cores=num.cores)
  return(Nodewise_Y_Xs_Xk_list)
}

### Sub-function: Treat one node as response 
### Return response vector, standardized design matrix and knockoff matrix
OneNode_func <- function(node_index, X, knockoff_method){
  
  X_design <- X[, -node_index]
  Y <- X[,node_index]
  
  X_Xk <- create.fixed(X_design, method=knockoff_method)   # function from "knockoff" package
  Xs <- X_Xk$X
  Xk <- X_Xk$Xk
  
  return(list(Y, Xs, Xk))
}
