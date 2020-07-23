### Compute the thresholds and obtain E_est (Algorithm 2 in the paper)

# input: 'W_matrix': feature statistic matri
#        'q': nominal FDR level
#        'offset': indicate FDR/mFDR control
#        'rule': AND or OR
#        '(a,c_1)': used in optimization problem
# output: 'E_est': estimated edge set

##########################################################
### Main function: 
E_est_givenW_func <- function(W_matrix, q, offset, rule, a){
  
  if(a==1) {c_a=1.93}
  if(a==0.01) {c_a=102}
  
  p <- dim(W_matrix)[1]
  W_sorted <- Sort_W_func(W_matrix)
  
  if(rule=="AND"){
    m_max <- floor( q*(p-1)/c_a - a*offset)
  }
  if(rule=="OR"){
    m_max <- floor( 1/2 * q*(p-1)/c_a - a*offset)
  }
  
  # if m_max < 0, will never find a feasible point, so thresholds=Inf
  if(m_max < 0){
    Finial_thre_vector <- rep(Inf,p)
  }
  
  # if m_max >= 0
  if(m_max >= 0){
    T_candidate_mat <- T_candidate_mat_func(W_sorted, m_max)
    
    ### Check the feasiability using each column of 'T_candidate_mat'
    Feasbile_indicate <- 0
    for (i in (m_max+1):1) {
      
      if(nrow(T_candidate_mat) ==1){    # R programming issue: for the case that 'T_candidate_mat' is a vector (i.e. m_max=0).
        thre_vector_temp <- as.vector(T_candidate_mat)
      }else{
        thre_vector_temp <- T_candidate_mat[,i]
      }
      
      constraint_vector_temp <- Constraint_vector_func(W_matrix, thre_vector_temp, rule, offset, a)
      
      # check if feasible or not
      if(rule=="OR") { 
        num.vio <- sum(constraint_vector_temp > q/(p*c_a))
      }
      if(rule=="AND") { 
        num.vio <- sum(constraint_vector_temp > 2*q/(p*c_a) )
      }
      if(num.vio==0) {Feasbile_indicate<-1; break}
    }
    
    if(Feasbile_indicate==0) {Finial_thre_vector <- rep(Inf,p)} else {Finial_thre_vector <- thre_vector_temp}
  }

  ### Recorve the E_est based on "W_matrix" and "Finial_thre_vector"
  E_temp <- apply(W_matrix, 2, `>=` , Finial_thre_vector)+0
  if(rule=="OR"){
    E_est <- ( (t(E_temp)+E_temp) != 0 )+0
  }
  if(rule=="AND"){
    E_est <- t(E_temp)*E_temp
  }
  
  return(E_est)
}

####################################################################################################################
####################################################################################################################
### Some functions used in the above main function:

### 1. Calculate the number of est-edges given threshold vector
E_num_func <- function(W_matrix, thre_vector, rule){
  
  E_selected <- apply(W_matrix, 2, `>=` , thre_vector)+0
  
  if(rule=="OR"){
    E_OR <- ( (t(E_selected)+E_selected) != 0 )+0
    E_num <- sum(E_OR)/2
  }
  if(rule=="AND"){
    E_AND <- t(E_selected)*E_selected
    E_num <- sum(E_AND)/2
  }
  
  return(E_num)
}

### 2. Calculate the constrains vector given threhsold vector 
Constraint_vector_func <- function(W_matrix, thre_vector, rule, offset, a){
  
  ### calculate demoninator |\hat(E)|
  E_num <- E_num_func(W_matrix, thre_vector, rule)
  
  constraint_vector <- mapply(constraint_func, 1:dim(W_matrix)[1], MoreArgs=list(W_matrix, thre_vector, offset, E_num, a) )
  
  return(constraint_vector)
}
# sub-func: calculate the constraint value for each row (in paper for each column)
constraint_func <- function(i, W_matrix, thre_vector, offset, E_num, a){
  V_neg_i <- sum(W_matrix[i,] <= -thre_vector[i])
  return( (a*offset + V_neg_i)/max(1,E_num) )
}

### 3. Return a W such that each row is in the abs-order but with signs (in paper each column)
Sort_W_func <- function(W_matrix){
  
  ### function for one row
  sort_onerow <- function(i, W_matrix){
    
    target_vector <- W_matrix[i,]
    target_vector_sorted <- sort(target_vector)                                      # sort in increasing for original values
    target_vector_abs_sorted <- sort(abs(target_vector_sorted), decreasing=TRUE)     # sort in decreasing for absolute values
    
    sign_location <- order(abs(target_vector_sorted), decreasing=TRUE)
    sign_sorted <- sign(target_vector_sorted)[sign_location]
    
    return(target_vector_abs_sorted*sign_sorted)        # return the vector in the abs order but with signs
  }
  
  W_sorted <- t(mapply(sort_onerow, 1:dim(W_matrix)[1], MoreArgs=list(W_matrix)))
  return(W_sorted)
}

### 4. Return a candidate matrix of thresholds for each row 
### T_candidate_mat: m_max+1 columns, p rows
### Each column relates to allowing 0, 1, 2, ..., m_max (m_max >= 0) negative values (i.e. |V_neg|)
T_candidate_mat_func <- function(W_sorted, m_max){
  
  T_candidate_mat <- t(mapply(One_row_T_candidate_func, 1:dim(W_sorted)[1], MoreArgs=list(W_sorted, m_max)))
  return(T_candidate_mat)
}

##### find T_candidate for one row of W
One_row_T_candidate_func <- function(i, W_sorted, m_max){
  
  target_vec <- W_sorted[i,]
  
  if(sum(target_vec>0)==0) {
    # if no positive value, return Inf for all m=0,1,...,m_max.
    return(rep(Inf, m_max+1))
  } else{
    # if there is at least one positive values and no negative values.
    # we just need to take the smallest positive value as threshold for all m=0,..,m_max.
    if (sum(target_vec<0)==0){
      return(rep(min(target_vec[which(target_vec>0)]), m_max+1))
    } 
  }
  
  # the case that there are at least one positive and one negative value
  
  T_candidate <- rep(Inf, m_max+1)  # each element corresponds to allowing 0, 1, 2, ..., m_max negative values
  
  # when allow 0 negtives (|V_neg|=0):
  if(target_vec[1]>0) {T_candidate[1] <- target_vec[which(target_vec<=0)[1]-1]} else {T_candidate[1] <- Inf}
  
  # if m_max=0, we return the above in which only one element in T_candidate
  # if m_max >= 1
  if(m_max >= 1){
    # when allow i (=1,2,...,m_max) negtives.
    # we store the threshold at T[i+1] as we already stored the case when allowing 0 negtives
    for (i in 1:m_max) {
      # if the allowed negtives i is larger than the total number of negatives, we just need to take the threshold as 
      # allowing i-1 negatives
      if(i > sum(target_vec<0)){   # note that from above, sum(target_vec<0) >= 1, so i=1 won't satisfies this 'if'
        T_candidate[i+1] <- T_candidate[i]
      }else{ 
        if(i == sum(target_vec<0)){     # in this case, the following elements would be non-negative
          rest_vec <- target_vec[(i+1):length(target_vec)]
          if(sum(rest_vec>0)>0){  # if there is at least one positive value in the rest
            T_candidate[i+1] <- min(rest_vec[rest_vec>0])
          } else{ # otherwise the rest value would be 0, then it's ok to take 'T_candidate[i]'
            T_candidate[i+1] <- T_candidate[i]
          }
          # here we omit the case where i=length(target_vec) becasue it is unlikely to happen
        }else{      # when i < sum(target_vec<0) = length(which(target_vec<0)), so there is another negative behind it.
          if(target_vec[which(target_vec<0)[i]+1] > 0){  # if the next element of the ith negative is positive
            rest_vec <- target_vec[(i+1):length(target_vec)]
            T_candidate[i+1] <- target_vec[which(target_vec<0)[i+1]-1]
          }else{ # if the next element of the ith negative is negative
            T_candidate[i+1] <- target_vec[which(target_vec<0)[i]]
          }
        }
      }
    } # end for
  }
  return(abs(T_candidate))
}