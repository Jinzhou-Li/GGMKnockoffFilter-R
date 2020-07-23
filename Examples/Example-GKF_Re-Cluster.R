library(mvtnorm)     # generate multivariate Gaussian date

source("OtherFuncs/Graphs.R")
source("OtherFuncs/FdrPowerGraphFunc.R")

source("GGMknockoffFilter/Nodewise_Y_Xs_Xk.R")
source("GGMknockoffFilter/Z_max_lambda.R")
source("GGMknockoffFilter/Z_coef.R")
source("GGMknockoffFilter/Z_stat_func_list_ElasticNet.R")
source("GGMknockoffFilter/E_est_givenW_func.R")
source("GGMknockoffFilter/GKF_Re.R")

#################################################
n = 3000          # number of observations
q = 0.2           # nominal FDR level

### Cluster graph
p_block <- 40           # number of variables in each cluster
num_block <- 5
prob = 0.5
range = c(0.2,0.6)
add_minEigen <- 0.5
signs <- "Random"         

Result <- Cluster_ER_graph(num_block, p_block, prob, range, add_minEigen, signs)
Omega <- Result[[1]]
Sigma <- Result[[2]]

mu = rep(0,p_block*num_block)

### Generate data
set.seed(666)
X <- rmvnorm(n, mu, Sigma)

### Estimate the edge set using GKF with sample-splitting-recycling to choose hyperparameter
E_est <- GKF_Re(X, q, num.cores=10)

print(Fdp_Power_Graph_func(E_est, Omega))
