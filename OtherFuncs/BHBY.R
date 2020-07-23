######## BH/BY combined with p-values obtained from testing zero partial correlations
library(ppcor)

Multi_Test_Graph_func <- function(X, q, method){
  
  p.value.matrix <- pcor(X)$p.value
  # p-values corresponding to each edge since p.value.matrix is symmetric
  p.value.vector <- p.value.matrix[which(upper.tri(p.value.matrix)==T)]    
  
  p.dim <- dim(p.value.matrix)[1]
  num.edge <- p.dim*(p.dim-1)/2          # number of tests/edges
  
  index.matrix <- matrix(1:(p.dim^2), nrow=p.dim, ncol=p.dim)
  index.vector <- index.matrix[which(upper.tri(index.matrix)==T)]
  
  test.result <- p.adjust(p.value.vector, method)<=q     # significant ones (TRUE <-> Rej)
  
  # transform the multiple testing result to estimated edge
  E_est <- matrix(0,p.dim,p.dim)
  
  for(i in 1:num.edge){
    row.index <- index.vector[i] %% p.dim
    col.index <- ceiling(index.vector[i]/p.dim)
    E_est[row.index, col.index] <- E_est[col.index, row.index] <- as.numeric(test.result[i])
  }
  
  return(E_est)
}

