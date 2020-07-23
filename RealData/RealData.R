############################################################################
###################   Load packages and functions     ######################
############################################################################
library(huge) # non-paranomral transformation
library(SILGGM)    # implement GFC-L and GFC-SL

source("OtherFuncs/BHBY.R")
source("OtherFuncs/KO2.R")

source("GGMknockoffFilter/Nodewise_Y_Xs_Xk.R")
source("GGMknockoffFilter/Z_max_lambda.R")
source("GGMknockoffFilter/Z_coef.R")
source("GGMknockoffFilter/Z_stat_func_list_ElasticNet.R")
source("GGMknockoffFilter/E_est_givenW_func.R")
source("GGMknockoffFilter/GKF_Re.R")

#################################################################################
# This public real data set can be downloaded from: https://journals.plos.org/ploscompbiol/article?id=%2010.1371/journal.pcbi.1006369 (see S2 Appendix).
load("RealData/sc_pan_T.RData")

##### select "num_var" variables with the largest variances
num_var <- 50
num.edge <- num_var*(num_var-1)/2

q <- 0.2

## Transformation: 1. log2 ; 2. nonparanormal: use huge.npn() function
sc_pan_T <- log2(sc_pan_T+1)
data <- huge.npn(sc_pan_T)
dim(data)

## 
var_vector <- apply(data, 2, var)
SelectedOnes <- order(var_vector, decreasing = TRUE)[1:num_var]
data_selected <- data[, SelectedOnes]

X <- data_selected

#########################################################################################################
########################## Knockoff based methods

Simulation_func <- function(index){
  print(index)
  
  ##### 1. Sample-splitting-recycling with GKF
  E_GKF_1 <- GKF_Re(X, q)

  return(list(E_GKF_1))
}

#########################################################################################################
Repli = 20
set.seed(999)

E_GKF_list <- t(mcmapply(Simulation_func, 1:Repli, mc.cores=20))

######################## other methods
E_est_BY <- Multi_Test_Graph_func(X, q, method="BY")
E_est_BH <- Multi_Test_Graph_func(X, q, method="BH")
E_est_L <- SILGGM(X, method="GFC_L", alpha=q)$global_decision[[1]]
E_est_SL <- SILGGM(X, method="GFC_SL", alpha=q)$global_decision[[1]]
ko.est2 <- GraphEstimation(X, q, plus = TRUE)
E_est_KO2 <- ko.est2$adjacency.matrix

Save_Result_list <- list(E_GKF_list, E_est_BY, E_est_BH, E_est_L, E_est_SL, E_est_KO2)

########################################################################################################
###################################### save data
save(Save_Result_list, file = "RealData/RealDataResults.RData")
