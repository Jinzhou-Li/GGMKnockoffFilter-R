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
p = 200
q = 0.2           # nominal FDR level

### Band graph
k= 10
c = 10
add_minEigen <- 0.5
b <- -0.6

Result <- Band_graph(p, k, b, c, add_minEigen)
Omega <- Result[[1]]
Sigma <- Result[[2]]

mu = rep(0,p)

### Generate data
set.seed(666)
X <- rmvnorm(n, mu, Sigma)

### Estimate the edge set using GKF with sample-splitting-recycling to choose hyperparameter
E_est <- GKF_Re(X, q, num.cores=10)

print(Fdp_Power_Graph_func(E_est, Omega))

