library(mvtnorm)     # generate multivariate Gaussian date

source("OtherFuncs/Graphs.R")
source("OtherFuncs/FdrPowerGraphFunc.R")

source("GGMknockoffFilter/Nodewise_Y_Xs_Xk.R")
source("GGMknockoffFilter/Z_max_lambda.R")
source("GGMknockoffFilter/Z_coef.R")
source("GGMknockoffFilter/E_est_givenW_func.R")
source("GGMknockoffFilter/GKF.R")

#################################################
n = 3000          # number of observations
p = 200           # number of variables
q = 0.2           # nominal FDR level

### Band graph
b <- -0.6
k= 10
c = 10
add_minEigen <- 0.5

Result <- Band_graph(p, k, b, c, add_minEigen)
Omega <- Result[[1]]
Sigma <- Result[[2]]

mu = rep(0,p)

### Generate data
set.seed(666)
X <- rmvnorm(n, mu, Sigma)

### Estimate the edge set using GKF with fixed hyperparameter
E_est <- GKF(X, q, offset=1, Z_matrix_function=Z_max_lambda, alpha=1,
             knockoff_method="equi", 
             rule="OR", a=0.01,
             W_matrix_BasedOnZ="max_sign", 
             num.cores=10)

print(Fdp_Power_Graph_func(E_est, Omega))

