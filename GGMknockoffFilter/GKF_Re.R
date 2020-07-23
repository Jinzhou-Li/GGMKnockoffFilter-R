### GGM knockoff filter with sample-splitting-recycling (Algprithm 4 in the paper)

# input: 'X': data matrix (n \times p)
#        'q': nominal FDR level
#        'offset': indicate FDR (offset=1) / mFDR (offset=0) control, the default value is 1
#        'knockoff_method_set': set of method to construct knockoffs
#        'rule_set': set of rule used to recover the estimated edge set from the estimated neighborhoods
#        'a_vec': set of parameter used in the optimization problem
#        'W_matrix_BasedOnZ_set': set of method used to construct feature statsitics (i.e. W-statistics) based on Z-statistics
#        'Z_stat_function_list': set of function used to construct Z-statistics (measure the importance of variable/knockoff to the response)
#        'num.cores': number of cores used for the parallel, the default value is 1
# output: 'E_est': estimated edge set

### The default hyperparameter space is set as in the simulation part of the paper. 

##########################################################

library(knockoff)        # for nodewisely constructing knockoffs
library(parallel)        # parallel when nodewisely construct knockoffs and Z-statistics

##########################################################
### Main function:
GKF_Re <- function(X, q, offset=1, Recycling=TRUE,
                   knockoff_method_set=c("equi", "sdp"), 
                   rule_set=c("AND","OR"), 
                   a_vec=c(1,0.01),
                   W_matrix_BasedOnZ_set=c("difference", "max_sign"), 
                   Z_stat_function_list = Z_stat_func_list_ElasticNet(alpha_vec=seq(0.2,1,0.2), lambda_quantile_vec=0.1*(1:10)),
                   num.cores=1){
  # split sample
  sam_index <- sample(1:nrow(X), round(nrow(X)/2), replace=FALSE)
  X_1 <- X[sam_index, ]
  X_2 <- X[-sam_index, ]
  
  #######################
  ### 1. Select hyperparameter based on X_1
  Hyperpara_trace_matrix <- NULL
  
  for (i_1 in 1:length(knockoff_method_set)) {
    Nodewise_Y_Xs_Xk_list <- Nodewise_Y_Xs_Xk(X_1, knockoff_method_set[i_1], num.cores)
    
    for (i_2 in 1:length(Z_stat_function_list)) {
      
      Z_matrix_function <- Z_stat_function_list[[i_2]]
      
      Z_return <- Z_matrix_function(Nodewise_Y_Xs_Xk_list, num.cores)
      Z_ori_return <- Z_return[[1]]
      Z_knock_return <- Z_return[[2]]
      
      # if 'Z_matrix_function' contains many functions and 'Z_ori_return' is an array
      if(length(dim(Z_ori_return))==3){
        for(i_2_sub in 1:dim(Z_ori_return)[3]){   # extra for-loop compared to below
          
          Z_ori_matrix <- Z_ori_return[,,i_2_sub]
          Z_knock_matrix <- Z_knock_return[,,i_2_sub]
          
          for (i_3 in 1:length(W_matrix_BasedOnZ_set)) {
            if(W_matrix_BasedOnZ_set[i_3]=="difference"){
              W_matrix = Z_ori_matrix - Z_knock_matrix
            }
            if(W_matrix_BasedOnZ_set[i_3]=="max_sign"){
              W_matrix = sign(Z_ori_matrix-Z_knock_matrix) * pmax(Z_ori_matrix,Z_knock_matrix)
            }
            
            if(sum(diag(W_matrix))!=0) stop("sum(diag(W_matrix))!=0 !")
            
            for (i_4 in 1:length(rule_set)) {
              for (i_5 in 1:length(a_vec)) {
                E_est <- E_est_givenW_func(W_matrix, q, offset, rule_set[i_4], a_vec[i_5])
                if(sum(diag(E_est))!=0) stop("sum(diag(E_est))!=0 !")
                
                num_disc <- sum(E_est)/2
                
                Hyperpara_trace_matrix <- rbind(Hyperpara_trace_matrix, c(i_1, i_2, i_2_sub, i_3, i_4, i_5, num_disc))
              } # end i_5
            }  # end i_4
          }  # end i_3
        } # end i_2_sub
      } # end if
      
      # if 'Z_matrix_function' contains only one function and 'Z_ori_return' is a matrix
      if(length(dim(Z_ori_return))==2){
        i_2_sub <- 0
        Z_ori_matrix <- Z_ori_return
        Z_knock_matrix <- Z_knock_return
        
        for (i_3 in 1:length(W_matrix_BasedOnZ_set)) {
          if(W_matrix_BasedOnZ_set[i_3]=="difference"){
            W_matrix = Z_ori_matrix - Z_knock_matrix
          }
          if(W_matrix_BasedOnZ_set[i_3]=="max_sign"){
            W_matrix = sign(Z_ori_matrix-Z_knock_matrix) * pmax(Z_ori_matrix,Z_knock_matrix)
          }
          if(sum(diag(W_matrix))!=0) stop("sum(diag(W_matrix))!=0 !")
          
          for (i_4 in 1:length(rule_set)) {
            for (i_5 in 1:length(a_vec)) {
              E_est <- E_est_givenW_func(W_matrix, q, offset, rule_set[i_4], a_vec[i_5])
              if(sum(diag(E_est))!=0) stop("sum(diag(E_est))!=0 !")
              
              num_disc <- sum(E_est)/2
              
              Hyperpara_trace_matrix <- rbind(Hyperpara_trace_matrix, c(i_1, i_2, i_2_sub, i_3, i_4, i_5, num_disc))
            } # end i_5
          }  # end i_4
        }  # end i_3
      } # end if
      
    }# end i_2
  } # end i_1
  
  # select the one with the maximal discoveries (randomly choose one if there is a tie)
  max_index_vec <- which(Hyperpara_trace_matrix[,7] == max(Hyperpara_trace_matrix[,7]))
  Hyperpara_selec_index <- max_index_vec[sample(length(max_index_vec), 1)]
  Hyperpara_selec <- Hyperpara_trace_matrix[Hyperpara_selec_index,]
  
  #######################
  ### 2. Implement the selected algorithm to X2 to get E_est
  
  Nodewise_Y_Xs_Xk_list2 <- Nodewise_Y_Xs_Xk(X_2, knockoff_method_set[Hyperpara_selec[1]], num.cores)
  Nodewise_Y_Xs_Xk_list_final <- Nodewise_Y_Xs_Xk_list2 
  
  # data recycling
  if(Recycling==TRUE){
    Nodewise_Y_list2 <- Nodewise_Y_Xs_Xk_list2[1,]
    Nodewise_Xs_list2 <- Nodewise_Y_Xs_Xk_list2[2,]
    Nodewise_Xk_list2 <- Nodewise_Y_Xs_Xk_list2[3,]
    
    # actually can delete the followoing line as 'Nodewise_Y_list1' and 'Nodewise_Y_list1' would be the same for any 'knockoff_method'
    Nodewise_Y_Xs_Xk_list1 <- Nodewise_Y_Xs_Xk(X_1, knockoff_method_set[Hyperpara_selec[1]], num.cores)  
    Nodewise_Y_list1 <- Nodewise_Y_Xs_Xk_list1[1,]
    Nodewise_Xs_list1 <- Nodewise_Y_Xs_Xk_list1[2,]
    
    Nodewise_Y_list_Re <- mapply(function(i, A_list, B_list){return(c(A_list[[i]], B_list[[i]]))}, 1:ncol(X), 
                                 MoreArgs=list(Nodewise_Y_list1, Nodewise_Y_list2), SIMPLIFY=FALSE)
    Nodewise_Xs_list_Re <- mapply(function(i, A_list, B_list){return(rbind(A_list[[i]], B_list[[i]]))}, 1:ncol(X),
                                  MoreArgs=list(Nodewise_Xs_list1, Nodewise_Xs_list2), SIMPLIFY=FALSE)
    Nodewise_Xk_list_Re <- mapply(function(i, A_list, B_list){return(rbind(A_list[[i]], B_list[[i]]))}, 1:ncol(X), 
                                  MoreArgs=list(Nodewise_Xs_list1, Nodewise_Xk_list2), SIMPLIFY=FALSE)
    Nodewise_Y_Xs_Xk_list_final[1,] <- Nodewise_Y_list_Re
    Nodewise_Y_Xs_Xk_list_final[2,] <- Nodewise_Xs_list_Re
    Nodewise_Y_Xs_Xk_list_final[3,] <- Nodewise_Xk_list_Re
  }
  
  Z_matrix_function <- Z_stat_function_list[[Hyperpara_selec[2]]]
  
  Z_return <- Z_matrix_function(Nodewise_Y_Xs_Xk_list_final, num.cores)
  Z_ori_return <- Z_return[[1]]
  Z_knock_return <- Z_return[[2]]
  
  if(Hyperpara_selec[3]==0){
    Z_ori_matrix <- Z_ori_return
    Z_knock_matrix <- Z_knock_return
  }
  if(Hyperpara_selec[3]!=0){
    Z_ori_matrix <- Z_ori_return[,,Hyperpara_selec[3]]
    Z_knock_matrix <- Z_knock_return[,,Hyperpara_selec[3]]
  }
  
  if(W_matrix_BasedOnZ_set[Hyperpara_selec[4]]=="difference"){
    W_matrix = Z_ori_matrix - Z_knock_matrix
  }
  if(W_matrix_BasedOnZ_set[Hyperpara_selec[4]]=="max_sign"){
    W_matrix = sign(Z_ori_matrix-Z_knock_matrix) * pmax(Z_ori_matrix,Z_knock_matrix)
  }
  
  E_est <- E_est_givenW_func(W_matrix, q, offset, rule_set[Hyperpara_selec[5]], a_vec[Hyperpara_selec[6]])
  
  return(E_est)
} # end function
