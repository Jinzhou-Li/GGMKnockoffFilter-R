######################################################################################
################ Band graph
Band_graph <- function(p, k, a, c, add_minEigen){
  
  Omega <- I <- diag(1,p)

  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(abs(i-j)<=k) Omega[i,j] <- Omega[j,i] <- sign(a)*( abs(a)^(abs((i-j)/c)) )
    }
  }
  
  permut <- sample(1:p, replace=FALSE)
  Omega <- as.matrix(Omega[permut,permut])

  # To make sure that Omega is positive semi-definite and the positive eigen-values are not too small
  lambda.min.abs <- abs(min(eigen(Omega)$values))
  Omega <- Omega + (lambda.min.abs + add_minEigen)*I
  Sigma <- solve(Omega)

  return( list(Omega, Sigma))
}

#######################################################################################
########### ER graph: entries sample from (-b,-a) (a,b)
ER_graph <- function(p, prob, range, add_minEigen, signs){
  
  Omega <- I <- diag(p)
  min <- min(range)
  max <- max(range)
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      u <- runif(1, min, max)
      delta <- rbinom(1, 1, prob)
      if(signs=="Positive"){Sign <- 1}
      if(signs=="Negative"){Sign <- -1}
      if(signs=="Random"){Sign <- sample(c(1,-1),1)}
      Omega[i,j] <- Omega[j,i] <- u*delta*Sign
    }
  }
  
  # To make sure that Omega is positive semi-definite and the positive eigen-values are not too small
  lambda.min.abs <- abs(min(eigen(Omega)$values));
  Omega <- Omega + (lambda.min.abs + add_minEigen)*I
  Sigma <- solve(Omega)
  
  return( list(Omega, Sigma))
}

#######################################################################################
Cluster_ER_graph <- function(num_cluster, p_cluster, prob, range, add_minEigen, signs){
  
  min <- min(range)
  max <- max(range)
  
  Omega <- diag(p_cluster)
  for(i in 1:(p_cluster-1)){
    for(j in (i+1):p_cluster){
      u <- runif(1, min, max)
      delta <- rbinom(1, 1, prob)
      if(signs=="Positive"){Sign <- 1}
      if(signs=="Negative"){Sign <- -1}
      if(signs=="Random"){Sign <- sample(c(1,-1),1)}
      Omega[i,j] <- Omega[j,i] <- u*delta*Sign
    }
  }
  
  if(num_cluster > 1){
    for (k in 2:num_cluster) {
      Omega_new <- diag(p_cluster)
      for(i in 1:(p_cluster-1)){
        for(j in (i+1):p_cluster){
          u <- runif(1, min, max)
          delta <- rbinom(1, 1, prob)
          if(signs=="Positive"){Sign <- 1}
          if(signs=="Negative"){Sign <- -1}
          if(signs=="Random"){Sign <- sample(c(1,-1),1)}
          Omega_new[i,j] <- Omega_new[j,i] <- u*delta*Sign
        }
      }
      Omega <- bdiag(Omega, Omega_new)
    }
  } # end if
  
  permut <- sample(1:(num_cluster*p_cluster), replace=FALSE)
  Omega <- as.matrix(Omega[permut,permut])
  
  # To make sure that Omega is positive semi-definite and the positive eigen-values are not too small
  I <- diag(p_cluster*num_cluster)
  lambda.min.abs <- abs(min(eigen(Omega)$values));
  Omega <- Omega + (lambda.min.abs + add_minEigen)*I
  Sigma <- solve(Omega)
  
  return( list(Omega, Sigma))
}

# Block_graph can be generated using 'Cluster_ER_graph' with e.g. prob=1, range=c(0.5,0.5)
